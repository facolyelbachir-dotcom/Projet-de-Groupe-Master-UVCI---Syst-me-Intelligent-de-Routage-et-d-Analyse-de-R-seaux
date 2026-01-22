
#include "backtracking.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* ==========================
 * Outils
 * ========================== */

static int in_set(const int* arr, size_t n, int v) {
    if (!arr) return 0;
    for (size_t i = 0; i < n; ++i) if (arr[i] == v) return 1;
    return 0;
}

/* Construit un graphe renversé (copie des latences/couts/bw/securite avec arcs inversés) */
static Graphe* graphe_renverse(const Graphe* g) {
    Graphe* r = graphe_creer(g->nb_noeuds);
    for (int u = 0; u < g->nb_noeuds; ++u) {
        for (Arete* a = g->noeuds[u].aretes; a; a = a->suivant) {
            graphe_ajouter_arete(r, a->destination, u, a->latence, a->bande_passante, a->cout, a->securite, 1);
        }
    }
    return r;
}

/* Dijkstra (latence) avec blocage de nœuds (block[v]!=0) et éventuellement d’arêtes non-sécurisées */
static void dijkstra_blocage(const Graphe* g, int src, float* out_dist,
                             const char* block_node, int interdire_non_secure) {
    int n = g->nb_noeuds;
    for (int i = 0; i < n; ++i) out_dist[i] = INFINITY;
    if (block_node && block_node[src]) return;

    int* prev = (int*)malloc(sizeof(int) * (size_t)n);
    char* vis = (char*)calloc((size_t)n, 1);
    for (int i = 0; i < n; ++i) prev[i] = -1;

    FileAttente* pq = file_creer(-1);
    out_dist[src] = 0.0f;
    file_push(pq, src, 0.0f);

    while (!file_est_vide(pq)) {
        int u = file_pop(pq);
        if (vis[u]) continue;
        vis[u] = 1;

        for (Arete* a = g->noeuds[u].aretes; a; a = a->suivant) {
            if (interdire_non_secure && a->securite == 0) continue;
            int v = a->destination;
            if (block_node && block_node[v]) continue;
            float alt = out_dist[u] + a->latence;
            if (alt < out_dist[v]) {
                out_dist[v] = alt; prev[v] = u;
                file_decrease_key(pq, v, alt);
            }
        }
    }

    file_detruire(pq);
    free(prev); free(vis);
}

/* ==========================
 *   Backtracking
 * ========================== */

typedef struct BTState {
    const Graphe* g;
    const Contraintes* cs;
    int src, dst;

    char* visite;                /* anti-cycles */
    int* chemin; size_t len; size_t cap;

    float latence;               /* cumul */
    float cout;                  /* cumul */
    float bw_min;                /* min cumulé */

    float* heur_lb;              /* Borne inférieure latence -> dst (par nœud) */

    Chemin best;                 /* meilleure solution */
} BTState;

static void bt_push(BTState* S, int v) {
    if (S->len == S->cap) {
        S->cap = S->cap ? S->cap * 2 : 16;
        S->chemin = (int*)realloc(S->chemin, sizeof(int) * S->cap);
    }
    S->chemin[S->len++] = v;
    S->visite[v] = 1;
}

static void bt_pop(BTState* S) {
    int v = S->chemin[S->len - 1];
    S->visite[v] = 0;
    S->len--;
}

static void bt_consider_best(BTState* S) {
    if (S->best.longueur == 0 || S->latence < S->best.latence_totale) {
        chemin_liberer(&S->best);
        S->best.noeuds = (int*)malloc(sizeof(int) * S->len);
        memcpy(S->best.noeuds, S->chemin, sizeof(int) * S->len);
        S->best.longueur = S->len;
        S->best.latence_totale = S->latence;
        S->best.cout_total = S->cout;
        S->best.bande_passante_min = S->bw_min;
    }
}

static int tous_must_visit_atteints(const BTState* S) {
    if (!S->cs->must_visit || S->cs->must_visit_count == 0) return 1;
    for (size_t i = 0; i < S->cs->must_visit_count; ++i) {
        int mv = S->cs->must_visit[i];
        int vu = 0;
        for (size_t t = 0; t < S->len; ++t) if (S->chemin[t] == mv) { vu = 1; break; }
        if (!vu) return 0;
    }
    return 1;
}

static void bt_dfs(BTState* S, int u) {
    /* Pruning sur coût / bande passante */
    if (S->cout > S->cs->cout_max) return;
    if (S->bw_min < S->cs->bande_passante_min) return;

    /* Heuristique LB (latence restante) */
    float lb = S->heur_lb[u];
    if (!isfinite(lb)) return; /* pas de chemin vers dst */
    if (S->best.longueur > 0 && (S->latence + lb) >= S->best.latence_totale) return;

    if (u == S->dst) {
        if (!tous_must_visit_atteints(S)) return;
        bt_consider_best(S);
        return;
    }

    /* Explorer voisins */
    for (Arete* a = S->g->noeuds[u].aretes; a; a = a->suivant) {
        if (S->cs->appliquer_securite && a->securite == 0) continue;
        int v = a->destination;
        if (S->visite[v]) continue;
        if (in_set(S->cs->exclude_nodes, S->cs->exclude_count, v)) continue;

        float old_lat = S->latence, old_c = S->cout, old_bw = S->bw_min;

        S->latence += a->latence;
        S->cout    += a->cout;
        if (S->len == 0) S->bw_min = a->bande_passante;
        else if (a->bande_passante < S->bw_min) S->bw_min = a->bande_passante;

        bt_push(S, v);
        bt_dfs(S, v);
        bt_pop(S);

        S->latence = old_lat; S->cout = old_c; S->bw_min = old_bw;
    }
}

/* prépare heuristique LB via Dijkstra sur graphe renversé depuis dst */
static void preparer_heuristique(const Graphe* g, int dst, const Contraintes* cs, float* out_lb) {
    Graphe* r = graphe_renverse(g);

    char* block = (char*)calloc((size_t)g->nb_noeuds, 1);
    for (size_t i = 0; i < cs->exclude_count; ++i) {
        int v = cs->exclude_nodes[i];
        if (v >= 0 && v < g->nb_noeuds) block[v] = 1;
    }
    dijkstra_blocage(r, dst, out_lb, block, cs->appliquer_securite);
    free(block);
    graphe_liberer(r);
}

Chemin chemin_contraint_backtracking(const Graphe* g, int src, int dst, const Contraintes* cs) {
    Chemin vide = chemin_vide();
    if (!g || !cs) return vide;
    if (src < 0 || src >= g->nb_noeuds || dst < 0 || dst >= g->nb_noeuds) return vide;
    if (in_set(cs->exclude_nodes, cs->exclude_count, src) ||
        in_set(cs->exclude_nodes, cs->exclude_count, dst)) return vide;

    BTState S;
    memset(&S, 0, sizeof(S));
    S.g = g; S.cs = cs; S.src = src; S.dst = dst;
    S.visite = (char*)calloc((size_t)g->nb_noeuds, 1);
    S.chemin = NULL; S.len = 0; S.cap = 0;
    S.latence = 0.0f; S.cout = 0.0f; S.bw_min = FLT_MAX;
    S.heur_lb = (float*)malloc(sizeof(float) * (size_t)g->nb_noeuds);
    S.best = chemin_vide();

    preparer_heuristique(g, dst, cs, S.heur_lb);

    /* point de départ */
    bt_push(&S, src);
    bt_dfs(&S, src);
    bt_pop(&S);

    free(S.visite); free(S.heur_lb);
    return S.best; /* peut être vide si aucune solution */
}
