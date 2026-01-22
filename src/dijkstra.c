
#include "dijkstra.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* ==========================
 * Helpers internes
 * ========================== */

static Arete* arete_nouvelle(int v, float lat, float bw, float cout, int sec) {
    Arete* a = (Arete*)malloc(sizeof(Arete));
    a->destination = v;
    a->latence = lat;
    a->bande_passante = bw;
    a->cout = cout;
    a->securite = sec ? 1 : 0;
    a->suivant = NULL;
    return a;
}

static float** matrice_creer(int n) {
    float** M = (float**)malloc(sizeof(float*) * (size_t)n);
    for (int i = 0; i < n; ++i) {
        M[i] = (float*)malloc(sizeof(float) * (size_t)n);
        for (int j = 0; j < n; ++j) M[i][j] = INFINITY;
        M[i][i] = 0.0f;
    }
    return M;
}

static void matrice_liberer(float** M, int n) {
    if (!M) return;
    for (int i = 0; i < n; ++i) free(M[i]);
    free(M);
}

/* ==========================
 *   API Graphe
 * ========================== */

Graphe* graphe_creer(int nb_noeuds) {
    Graphe* g = (Graphe*)malloc(sizeof(Graphe));
    g->nb_noeuds = nb_noeuds;
    g->noeuds = (Noeud*)calloc((size_t)nb_noeuds, sizeof(Noeud));
    for (int i = 0; i < nb_noeuds; ++i) {
        g->noeuds[i].id = i;
        g->noeuds[i].nom[0] = '\0';
        g->noeuds[i].aretes = NULL;
    }
    g->matrice_adjacence = matrice_creer(nb_noeuds);
    return g;
}

void graphe_liberer(Graphe* g) {
    if (!g) return;
    if (g->noeuds) {
        for (int i = 0; i < g->nb_noeuds; ++i) {
            Arete* a = g->noeuds[i].aretes;
            while (a) { Arete* nxt = a->suivant; free(a); a = nxt; }
        }
        free(g->noeuds);
    }
    matrice_liberer(g->matrice_adjacence, g->nb_noeuds);
    free(g);
}

void graphe_ajouter_arete(Graphe* g, int u, int v,
                          float latence, float bande_passante, float cout, int securite,
                          int directed) {
    if (!g || u < 0 || v < 0 || u >= g->nb_noeuds || v >= g->nb_noeuds) return;

    /* liste d’adjacence */
    Arete* a = arete_nouvelle(v, latence, bande_passante, cout, securite);
    a->suivant = g->noeuds[u].aretes;
    g->noeuds[u].aretes = a;

    /* matrice (latence) */
    if (latence < g->matrice_adjacence[u][v]) g->matrice_adjacence[u][v] = latence;

    if (!directed) {
        Arete* b = arete_nouvelle(u, latence, bande_passante, cout, securite);
        b->suivant = g->noeuds[v].aretes;
        g->noeuds[v].aretes = b;
        if (latence < g->matrice_adjacence[v][u]) g->matrice_adjacence[v][u] = latence;
    }
}

void graphe_supprimer_arete(Graphe* g, int u, int v) {
    if (!g || u < 0 || v < 0 || u >= g->nb_noeuds || v >= g->nb_noeuds) return;

    Arete* prev = NULL; Arete* cur = g->noeuds[u].aretes;
    while (cur) {
        if (cur->destination == v) {
            if (prev) prev->suivant = cur->suivant;
            else g->noeuds[u].aretes = cur->suivant;
            free(cur);
            break;
        }
        prev = cur; cur = cur->suivant;
    }

    /* remettre la cellule à +INF puis recalculer min latence si besoin */
    g->matrice_adjacence[u][v] = INFINITY;
    for (Arete* a2 = g->noeuds[u].aretes; a2; a2 = a2->suivant) {
        if (a2->destination == v && a2->latence < g->matrice_adjacence[u][v]) {
            g->matrice_adjacence[u][v] = a2->latence;
        }
    }
}

void graphe_set_nom(Graphe* g, int id, const char* nom) {
    if (!g || id < 0 || id >= g->nb_noeuds || !nom) return;
    strncpy(g->noeuds[id].nom, nom, sizeof(g->noeuds[id].nom)-1);
    g->noeuds[id].nom[sizeof(g->noeuds[id].nom)-1] = '\0';
}

void graphe_sync_matrice(Graphe* g) {
    if (!g) return;
    for (int i = 0; i < g->nb_noeuds; ++i) {
        for (int j = 0; j < g->nb_noeuds; ++j) {
            g->matrice_adjacence[i][j] = (i==j) ? 0.0f : INFINITY;
        }
    }
    for (int u = 0; u < g->nb_noeuds; ++u) {
        for (Arete* a = g->noeuds[u].aretes; a; a = a->suivant) {
            if (a->latence < g->matrice_adjacence[u][a->destination]) {
                g->matrice_adjacence[u][a->destination] = a->latence;
            }
        }
    }
}

/* ==========================
 *   API File de priorité
 * ========================== */

static Paquet* paquet_new(int id, float key) {
    Paquet* p = (Paquet*)malloc(sizeof(Paquet));
    p->id = id;
    p->taille_Mo = key;              /* on stocke la priorité flottante ici */
    p->priorite = (int)llroundf(key * 1000.0f); /* copie int (informative) */
    p->source = 0; p->destination = 0;
    p->precedent = NULL; p->suivant = NULL;
    return p;
}

FileAttente* file_creer(int capacite_max) {
    FileAttente* f = (FileAttente*)malloc(sizeof(FileAttente));
    f->tete = f->queue = NULL;
    f->taille_actuelle = 0;
    f->capacite_max = capacite_max;
    return f;
}

void file_detruire(FileAttente* f) {
    if (!f) return;
    Paquet* c = f->tete;
    while (c) { Paquet* n = c->suivant; free(c); c = n; }
    free(f);
}

int file_est_vide(const FileAttente* f) { return !f || f->tete == NULL; }

/* Insertion triée (priorité croissante) */
void file_push(FileAttente* f, int id, float key) {
    if (!f) return;
    if (f->capacite_max >= 0 && f->taille_actuelle >= f->capacite_max) return;

    Paquet* p = paquet_new(id, key);

    if (!f->tete) {
        f->tete = f->queue = p;
    } else if (key < f->tete->taille_Mo) {
        /* insertion en tête */
        p->suivant = f->tete;
        f->tete->precedent = p;
        f->tete = p;
    } else {
        /* insertion ordonnée */
        Paquet* cur = f->tete;
        while (cur->suivant && cur->suivant->taille_Mo <= key) cur = cur->suivant;
        p->suivant = cur->suivant;
        p->precedent = cur;
        if (cur->suivant) cur->suivant->precedent = p;
        cur->suivant = p;
        if (f->queue == cur) f->queue = p;
    }
    f->taille_actuelle++;
}

/* Retire la tête */
int file_pop(FileAttente* f) {
    if (!f || !f->tete) return -1;
    Paquet* p = f->tete;
    int id = p->id;
    f->tete = p->suivant;
    if (f->tete) f->tete->precedent = NULL; else f->queue = NULL;
    free(p);
    f->taille_actuelle--;
    return id;
}

/* Retire la première occurrence d'id s'il existe. */
static int file_remove_id(FileAttente* f, int id) {
    if (!f) return 0;
    Paquet* c = f->tete;
    while (c) {
        if (c->id == id) {
            if (c->precedent) c->precedent->suivant = c->suivant;
            else f->tete = c->suivant;
            if (c->suivant) c->suivant->precedent = c->precedent;
            else f->queue = c->precedent;
            free(c);
            f->taille_actuelle--;
            return 1;
        }
        c = c->suivant;
    }
    return 0;
}

void file_decrease_key(FileAttente* f, int id, float new_key) {
    if (!f) return;
    /* supprime si présent puis réinsère */
    file_remove_id(f, id);
    file_push(f, id, new_key);
}

/* ==========================
 *   Outils Chemin
 * ========================== */

Chemin chemin_vide(void) {
    Chemin c; c.noeuds = NULL; c.longueur = 0;
    c.latence_totale = 0.0f; c.cout_total = 0.0f; c.bande_passante_min = FLT_MAX;
    return c;
}

void chemin_liberer(Chemin* c) {
    if (!c) return;
    free(c->noeuds);
    c->noeuds = NULL;
    c->longueur = 0;
}

/* Construit un chemin (avec métriques) depuis prev[] */
static Chemin chemin_depuis_prev(const Graphe* g, int src, int dst, const int* prev, const float* dist) {
    Chemin C = chemin_vide();
    if (!g || !prev || !dist) return C;
    if (dst < 0 || dst >= g->nb_noeuds) return C;
    if (!isfinite(dist[dst])) return C;

    int cap = 16, len = 0;
    int* rev = (int*)malloc(sizeof(int) * (size_t)cap);
    for (int v = dst; v != -1; v = prev[v]) {
        if (len >= cap) { cap *= 2; rev = (int*)realloc(rev, sizeof(int) * (size_t)cap); }
        rev[len++] = v;
        if (v == src) break;
    }
    if (len == 0 || rev[len-1] != src) { free(rev); return C; }

    C.longueur = (size_t)len;
    C.noeuds = (int*)malloc(sizeof(int) * C.longueur);
    for (size_t i = 0; i < C.longueur; ++i) C.noeuds[i] = rev[len - 1 - (int)i];

    /* Calcul des métriques (latence connue) */
    C.latence_totale = dist[dst];
    C.cout_total = 0.0f;
    C.bande_passante_min = FLT_MAX;

    for (size_t i = 0; i + 1 < C.longueur; ++i) {
        int u = C.noeuds[i], v = C.noeuds[i+1];
        for (Arete* a = g->noeuds[u].aretes; a; a = a->suivant) {
            if (a->destination == v) {
                C.cout_total += a->cout;
                if (a->bande_passante < C.bande_passante_min)
                    C.bande_passante_min = a->bande_passante;
                break;
            }
        }
    }
    if (C.longueur <= 1) C.bande_passante_min = FLT_MAX;

    free(rev);
    return C;
}

/* ==========================
 *   Dijkstra
 * ========================== */

Chemin dijkstra_chemin(const Graphe* g, int src, int dst) {
    Chemin vide = chemin_vide();
    if (!g || src < 0 || dst < 0 || src >= g->nb_noeuds || dst >= g->nb_noeuds) return vide;

    float* dist = (float*)malloc(sizeof(float) * (size_t)g->nb_noeuds);
    int* prev   = (int*)malloc(sizeof(int) * (size_t)g->nb_noeuds);
    char* vis   = (char*)calloc((size_t)g->nb_noeuds, 1);

    for (int i = 0; i < g->nb_noeuds; ++i) { dist[i] = INFINITY; prev[i] = -1; }
    dist[src] = 0.0f;

    FileAttente* pq = file_creer(-1);
    file_push(pq, src, 0.0f);

    while (!file_est_vide(pq)) {
        int u = file_pop(pq);
        if (vis[u]) continue;
        vis[u] = 1;
        if (u == dst) break;

        for (Arete* a = g->noeuds[u].aretes; a; a = a->suivant) {
            int v = a->destination;
            float alt = dist[u] + a->latence;
            if (alt < dist[v]) {
                dist[v] = alt; prev[v] = u;
                file_decrease_key(pq, v, alt);
            }
        }
    }

    Chemin C = chemin_depuis_prev(g, src, dst, prev, dist);

    file_detruire(pq);
    free(dist); free(prev); free(vis);
    return C;
}

/* ==========================
 *   Bellman-Ford
 * ========================== */

Chemin bellman_ford_chemin(const Graphe* g, int src, int dst, int* has_negative_cycle) {
    Chemin vide = chemin_vide();
    if (has_negative_cycle) *has_negative_cycle = 0;
    if (!g || src < 0 || dst < 0 || src >= g->nb_noeuds || dst >= g->nb_noeuds) return vide;

    float* dist = (float*)malloc(sizeof(float) * (size_t)g->nb_noeuds);
    int* prev   = (int*)malloc(sizeof(int) * (size_t)g->nb_noeuds);

    for (int i = 0; i < g->nb_noeuds; ++i) { dist[i] = INFINITY; prev[i] = -1; }
    dist[src] = 0.0f;

    /* Relaxation V-1 fois */
    for (int it = 0; it < g->nb_noeuds - 1; ++it) {
        int changed = 0;
        for (int u = 0; u < g->nb_noeuds; ++u) {
            if (!isfinite(dist[u])) continue;
            for (Arete* a = g->noeuds[u].aretes; a; a = a->suivant) {
                int v = a->destination;
                float alt = dist[u] + a->latence;
                if (alt < dist[v]) {
                    dist[v] = alt; prev[v] = u; changed = 1;
                }
            }
        }
        if (!changed) break;
    }

    /* Détection cycle négatif atteignable */
    for (int u = 0; u < g->nb_noeuds; ++u) {
        if (!isfinite(dist[u])) continue;
        for (Arete* a = g->noeuds[u].aretes; a; a = a->suivant) {
            int v = a->destination;
            float alt = dist[u] + a->latence;
            if (alt < dist[v]) {
                if (has_negative_cycle) *has_negative_cycle = 1;
                free(dist); free(prev);
                return vide;
            }
        }
    }

    Chemin C = chemin_depuis_prev(g, src, dst, prev, dist);
    free(dist); free(prev);
    return C;
}

/* ==========================
 *   K plus courts chemins (Yen simplifié)
 * ========================== */

typedef struct CheminArray {
    Chemin* data;
    size_t  size;
    size_t  cap;
} CheminArray;

static void ca_init(CheminArray* A) { A->data = NULL; A->size = 0; A->cap = 0; }
static void ca_push(CheminArray* A, Chemin c) {
    if (A->size == A->cap) {
        A->cap = A->cap ? A->cap * 2 : 8;
        A->data = (Chemin*)realloc(A->data, sizeof(Chemin) * A->cap);
    }
    A->data[A->size++] = c;
}

/* clone superficiel du graphe (copie des listes) */
static Graphe* graphe_cloner(const Graphe* g) {
    Graphe* h = graphe_creer(g->nb_noeuds);
    for (int u = 0; u < g->nb_noeuds; ++u) {
        for (Arete* a = g->noeuds[u].aretes; a; a = a->suivant) {
            graphe_ajouter_arete(h, u, a->destination, a->latence, a->bande_passante, a->cout, a->securite, 1);
        }
    }
    return h;
}

static int chemins_compare_latence(const Chemin* a, const Chemin* b) {
    if (a->longueur == 0 && b->longueur > 0) return 0;
    if (b->longueur == 0 && a->longueur > 0) return 1;
    return a->latence_totale < b->latence_totale;
}

Chemin* k_plus_courts_chemins(const Graphe* g, int src, int dst, size_t K, size_t* out_count) {
    if (out_count) *out_count = 0;
    CheminArray A; ca_init(&A);
    CheminArray B; ca_init(&B);

    Chemin p0 = dijkstra_chemin(g, src, dst);
    if (p0.longueur == 0) {
        if (out_count) *out_count = 0;
        return NULL;
    }
    ca_push(&A, p0);

    for (size_t k = 1; k < K; ++k) {
        Chemin dernier = A.data[A.size - 1];

        for (size_t i = 0; i + 1 < dernier.longueur; ++i) {
            int spur_node = dernier.noeuds[i];

            /* clone et suppression des arêtes reconstituant les mêmes préfixes */
            Graphe* temp = graphe_cloner(g);

            for (size_t a = 0; a < A.size; ++a) {
                Chemin* P = &A.data[a];
                if (P->longueur > i) {
                    int meme_prefixe = 1;
                    for (size_t t = 0; t <= i; ++t) {
                        if (P->noeuds[t] != dernier.noeuds[t]) { meme_prefixe = 0; break; }
                    }
                    if (meme_prefixe) {
                        int u = P->noeuds[i];
                        int v = P->noeuds[i+1];
                        graphe_supprimer_arete(temp, u, v);
                    }
                }
            }

            /* calcul du spur path depuis spur_node */
            Chemin spur = dijkstra_chemin(temp, spur_node, dst);
            if (spur.longueur > 0) {
                /* concat root[0..i] + spur[1..] */
                size_t root_len = i + 1;
                size_t new_len = root_len + (spur.longueur - 1);
                Chemin cand = chemin_vide();
                cand.noeuds = (int*)malloc(sizeof(int) * new_len);
                cand.longueur = new_len;

                for (size_t t = 0; t < root_len; ++t) cand.noeuds[t] = dernier.noeuds[t];
                for (size_t t = 1; t < spur.longueur; ++t) cand.noeuds[root_len + t - 1] = spur.noeuds[t];

                /* recalculer métriques sur g */
                cand.latence_totale = 0.0f;
                cand.cout_total = 0.0f;
                cand.bande_passante_min = FLT_MAX;
                int valid = 1;
                for (size_t t = 0; t + 1 < cand.longueur; ++t) {
                    int u = cand.noeuds[t], v = cand.noeuds[t+1];
                    int found = 0;
                    for (Arete* a2 = g->noeuds[u].aretes; a2; a2 = a2->suivant) {
                        if (a2->destination == v) {
                            found = 1;
                            cand.latence_totale += a2->latence;
                            cand.cout_total += a2->cout;
                            if (a2->bande_passante < cand.bande_passante_min)
                                cand.bande_passante_min = a2->bande_passante;
                            break;
                        }
                    }
                    if (!found) { valid = 0; break; }
                }
                if (valid) ca_push(&B, cand);
                else chemin_liberer(&cand);
            }
            chemin_liberer(&spur);
            graphe_liberer(temp);
        }

        if (B.size == 0) break;

        /* sélectionner le meilleur candidat */
        size_t best_idx = 0;
        for (size_t i = 1; i < B.size; ++i) {
            if (chemins_compare_latence(&B.data[i], &B.data[best_idx])) best_idx = i;
        }
        Chemin best = B.data[best_idx];
        /* remove best (swap-pop) */
        B.data[best_idx] = B.data[B.size - 1];
        B.size--;

        ca_push(&A, best);

        /* nettoyer le reste des candidats */
        for (size_t i = 0; i < B.size; ++i) chemin_liberer(&B.data[i]);
        free(B.data); ca_init(&B);
    }

    if (out_count) *out_count = A.size;
    return A.data; /* à libérer par l'appelant */
}

/* ==========================
 *   Indice de redondance
 * ========================== */

float indice_redondance(const Chemin* chemins, size_t nb, int consider_edges) {
    if (!chemins || nb <= 1) return 0.0f;
    double acc = 0.0;
    size_t pairs = 0;

    for (size_t i = 0; i < nb; ++i) {
        for (size_t j = i + 1; j < nb; ++j) {
            double inter = 0.0, base = 0.0;

            if (consider_edges) {
                int e1 = (int)(chemins[i].longueur > 0 ? (int)chemins[i].longueur - 1 : 0);
                int e2 = (int)(chemins[j].longueur > 0 ? (int)chemins[j].longueur - 1 : 0);
                base = (double)(e1 + e2);
                if (base <= 0.0) continue;
                for (size_t a = 0; a + 1 < chemins[i].longueur; ++a) {
                    int u1 = chemins[i].noeuds[a], v1 = chemins[i].noeuds[a+1];
                    for (size_t b = 0; b + 1 < chemins[j].longueur; ++b) {
                        int u2 = chemins[j].noeuds[b], v2 = chemins[j].noeuds[b+1];
                        if (u1 == u2 && v1 == v2) inter += 2.0;
                    }
                }
            } else {
                base = (double)(chemins[i].longueur + chemins[j].longueur);
                if (base <= 0.0) continue;
                for (size_t a = 0; a < chemins[i].longueur; ++a) {
                    int u = chemins[i].noeuds[a];
                    for (size_t b = 0; b < chemins[j].longueur; ++b) {
                        if (u == chemins[j].noeuds[b]) { inter += 2.0; break; }
                    }
                }
            }

            acc += inter / base;
            pairs++;
        }
    }
    return (pairs == 0) ? 0.0f : (float)(acc / (double)pairs);
}
