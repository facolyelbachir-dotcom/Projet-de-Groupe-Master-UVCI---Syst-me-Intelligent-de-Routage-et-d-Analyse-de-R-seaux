
#include "dijkstra.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>

#define INFINI_FLT (1e30f)

/* =========================
 *  File de priorité (liste chaînée triée par priorité croissante)
 * ========================= */

typedef struct PQNode {
    int   sommet;
    float priorite;      /* distance actuelle */
    struct PQNode* prev;
    struct PQNode* next;
} PQNode;

typedef struct {
    PQNode* head;        /* plus petite priorité en tête */
} PQueue;

static PQueue* pq_creer(void) {
    PQueue* q = (PQueue*)malloc(sizeof(PQueue));
    if (q) q->head = NULL;
    return q;
}

static void pq_detruire(PQueue* q) {
    if (!q) return;
    PQNode* cur = q->head;
    while (cur) {
        PQNode* n = cur->next;
        free(cur);
        cur = n;
    }
    free(q);
}

/* insertion triée O(n), pop en O(1) */
static void pq_push(PQueue* q, int sommet, float priorite) {
    PQNode* node = (PQNode*)malloc(sizeof(PQNode));
    if (!node) return; /* hors mémoire -> on laisse l'algorithme échouer plus loin */
    node->sommet = sommet;
    node->priorite = priorite;
    node->prev = node->next = NULL;

    if (!q->head) { q->head = node; return; }

    PQNode* cur = q->head;
    PQNode* prev = NULL;
    while (cur && cur->priorite <= priorite) {
        prev = cur;
        cur = cur->next;
    }
    if (!prev) { /* insertion en tête */
        node->next = q->head;
        q->head->prev = node;
        q->head = node;
    } else {
        node->next = prev->next;
        node->prev = prev;
        prev->next = node;
        if (node->next) node->next->prev = node;
    }
}

static bool pq_pop_min(PQueue* q, int* sommet, float* priorite) {
    if (!q || !q->head) return false;
    PQNode* node = q->head;
    q->head = node->next;
    if (q->head) q->head->prev = NULL;
    if (sommet)   *sommet = node->sommet;
    if (priorite) *priorite = node->priorite;
    free(node);
    return true;
}

static bool pq_est_vide(const PQueue* q) {
    return !q || !q->head;
}

/* =========================
 *  Helpers Graphe
 * ========================= */

static const Arete* trouver_arete(const Graphe* g, int u, int v) {
    if (!g || u < 0 || v < 0 || u >= g->nb_noeuds || v >= g->nb_noeuds) return NULL;
    const Arete* a = g->noeuds[u].aretes;
    while (a) {
        if (a->destination == v) return a;
        a = a->suivant;
    }
    return NULL;
}

static void calculer_mesures_chemin(const Graphe* g, const int* parent, int source, int destination,
                                    float* out_lat, float* out_cout, float* out_bpmin) {
    float lat = 0.0f, cout = 0.0f;
    float bpmin = INFINI_FLT;

    int v = destination;
    while (v != source && v != -1) {
        int u = parent[v];
        if (u == -1) break;
        const Arete* e = trouver_arete(g, u, v);
        if (e) {
            lat  += e->latence;
            cout += e->cout;
            if (e->bande_passante < bpmin) bpmin = e->bande_passante;
        }
        v = u;
    }

    if (out_lat)   *out_lat = lat;
    if (out_cout)  *out_cout = cout;
    if (out_bpmin) *out_bpmin = (bpmin == INFINI_FLT ? 0.0f : bpmin);
}

static Chemin* reconstruire_chemin(const Graphe* g, const int* parent, int source, int destination) {
    /* Compter la longueur */
    int len = 0;
    for (int v = destination; v != -1; v = parent[v]) {
        ++len;
        if (v == source) break;
    }
    if (len == 0) return NULL;

    Chemin* c = (Chemin*)malloc(sizeof(Chemin));
    if (!c) return NULL;
    c->noeuds = (int*)malloc(sizeof(int) * len);
    if (!c->noeuds) { free(c); return NULL; }
    c->longueur = len;
    c->latence_totale = 0.0f;
    c->cout_total = 0.0f;
    c->bande_passante_min = 0.0f;

    /* Remplir à l'envers puis renverser */
    int idx = len - 1;
    for (int v = destination; v != -1; v = parent[v]) {
        c->noeuds[idx--] = v;
        if (v == source) break;
    }

    /* Calcul des mesures */
    calculer_mesures_chemin(g, parent, source, destination,
                            &c->latence_totale, &c->cout_total, &c->bande_passante_min);
    return c;
}

/* =========================
 *  API publique
 * ========================= */

void chemin_liberer(Chemin* c) {
    if (!c) return;
    free(c->noeuds);
    free(c);
}

int dijkstra_remplir_tableaux(const Graphe* g, int source, float* dist, int* parent) {
    if (!g || !dist || !parent) return 1;
    if (source < 0 || source >= g->nb_noeuds) return 2;

    const int V = g->nb_noeuds;

    for (int i = 0; i < V; ++i) {
        dist[i] = INFINI_FLT;
        parent[i] = -1;
    }
    dist[source] = 0.0f;

    /* Vérification d'arêtes négatives (non supportées) */
    for (int u = 0; u < V; ++u) {
        for (const Arete* e = g->noeuds[u].aretes; e; e = e->suivant) {
            if (e->latence < 0.0f) {
                return 3; /* arête négative détectée */
            }
        }
    }

    PQueue* pq = pq_creer();
    if (!pq) return 4;
    pq_push(pq, source, 0.0f);

    while (!pq_est_vide(pq)) {
        int u; float du;
        pq_pop_min(pq, &u, &du);

        /* Entrée obsolète ? */
        if (du > dist[u]) continue;

        /* Relâchements */
        for (const Arete* e = g->noeuds[u].aretes; e; e = e->suivant) {
            int v = e->destination;
            float w = e->latence;
            /* Dijkstra exige w >= 0, déjà vérifié avant */
            float nv = dist[u] + w;
            if (nv < dist[v]) {
                dist[v] = nv;
                parent[v] = u;
                pq_push(pq, v, nv); /* pas de decrease-key: duplication + test d'obsolescence */
            }
        }
    }

    pq_detruire(pq);
    return 0;
}

Chemin* dijkstra_plus_court_par_latence(const Graphe* g, int source, int destination, int* code_retour) {
    if (code_retour) *code_retour = 0;
    if (!g) { if (code_retour) *code_retour = 1; return NULL; }
    if (source < 0 || destination < 0 || source >= g->nb_noeuds || destination >= g->nb_noeuds) {
        if (code_retour) *code_retour = 1;
        return NULL;
    }

    const int V = g->nb_noeuds;
    float* dist = (float*)malloc(sizeof(float) * V);
    int* parent = (int*)malloc(sizeof(int) * V);
    if (!dist || !parent) {
        free(dist); free(parent);
        if (code_retour) *code_retour = 4;
        return NULL;
    }

    int rc = dijkstra_remplir_tableaux(g, source, dist, parent);
    if (rc != 0) {
        if (code_retour) *code_retour = (rc == 3 ? 3 : 4);
        free(dist); free(parent);
        return NULL;
    }

    if (dist[destination] >= INFINI_FLT / 2) {
        if (code_retour) *code_retour = 2; /* pas de chemin */
        free(dist); free(parent);
        return NULL;
    }

    Chemin* c = reconstruire_chemin(g, parent, source, destination);
    free(dist);
    free(parent);
    return c;
}
``
