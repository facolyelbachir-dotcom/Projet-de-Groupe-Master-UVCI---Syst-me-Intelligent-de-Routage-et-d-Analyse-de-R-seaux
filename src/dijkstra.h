
#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <stddef.h>

/* ============================================================
 *   Représentations du graphe (OBLIGATOIRE : 2 représentations)
 * ============================================================ */

typedef struct Arete {
    int destination;             /* id du voisin */
    float latence;               /* poids principal (pour Dijkstra/BF) */
    float bande_passante;        /* contrainte (min sur chemin) */
    float cout;                  /* contrainte (somme sur chemin) */
    int   securite;              /* 1=secure, 0=non-secure (liaison) */
    struct Arete* suivant;       /* liste d’adjacence */
} Arete;

typedef struct Noeud {
    int id;
    char nom[50];
    Arete* aretes;               /* tête de la liste d’adjacence */
} Noeud;

typedef struct Graphe {
    int nb_noeuds;
    Noeud* noeuds;

    /* Représentation alternative : matrice d'adjacence (latences).
       Si pas d’arête i->j : matrice_adjacence[i][j] = +INF */
    float** matrice_adjacence;
} Graphe;

/* ==========================
 * File de priorité (OBLIG.)
 * ==========================
 * Implémentée par une liste chaînée double ordonnée par priorité croissante.
 * On utilise 'taille_Mo' pour stocker la priorité flottante (distance).
 * 'priorite' est conservé (int) pour satisfaire la contrainte de structure.
 */

typedef struct Paquet {
    int id;                 /* identifiant (ici : id de nœud) */
    int priorite;           /* copie entière optionnelle de la priorité */
    float taille_Mo;        /* priorité réelle (distance, plus petit = mieux) */
    int source;             /* inutilisé ici mais conservé */
    int destination;        /* inutilisé ici mais conservé */
    struct Paquet* precedent;
    struct Paquet* suivant;
} Paquet;

typedef struct FileAttente {
    Paquet* tete;
    Paquet* queue;
    int taille_actuelle;
    int capacite_max;       /* non restrictif ici (peut valoir -1 pour illimité) */
} FileAttente;

/* ==========================
 *   Chemin utilitaire
 * ========================== */
typedef struct Chemin {
    int* noeuds;                /* séquence de nœuds */
    size_t longueur;            /* nombre de nœuds */
    float latence_totale;       /* somme latences */
    float cout_total;           /* somme coûts */
    float bande_passante_min;   /* min bande passante le long du chemin */
} Chemin;

/* ==========================
 *      API Graphe
 * ========================== */

Graphe* graphe_creer(int nb_noeuds);
void    graphe_liberer(Graphe* g);

/* Ajoute arête u->v (directed=1), ou u<->v (directed=0). Met à jour la matrice (latence). */
void graphe_ajouter_arete(Graphe* g, int u, int v,
                          float latence, float bande_passante, float cout, int securite,
                          int directed);

/* Supprime arête u->v (si présente). */
void graphe_supprimer_arete(Graphe* g, int u, int v);

/* Renomme un nœud (optionnel, utile pour debug/affichage). */
void graphe_set_nom(Graphe* g, int id, const char* nom);

/* Réinitialise/Reconstruit la matrice d’adjacence (latences) depuis les listes. */
void graphe_sync_matrice(Graphe* g);

/* ==========================
 *   API File de priorité
 * ========================== */

FileAttente* file_creer(int capacite_max /* -1 = illimité */);
void         file_detruire(FileAttente* f);
int          file_est_vide(const FileAttente* f);

/* Insertion ordonnée (priorité = key). Complexité O(n). */
void         file_push(FileAttente* f, int id, float key);

/* Extrait et retourne l'id en tête (plus petite priorité). Retourne -1 si vide. */
int          file_pop(FileAttente* f);

/* Diminue la clé d’un élément (id). Si absent, insère. Complexité O(n). */
void         file_decrease_key(FileAttente* f, int id, float new_key);

/* ==========================
 *   Outils Chemin
 * ========================== */
Chemin chemin_vide(void);
void   chemin_liberer(Chemin* c);

/* ==========================
 *  Algorithmes de base
 * ========================== */

/* Dijkstra (latence comme poids) — file de priorité par liste double. */
Chemin dijkstra_chemin(const Graphe* g, int src, int dst);

/* Bellman-Ford (latence comme poids, gère poids négatifs).
   has_negative_cycle != NULL => *has_negative_cycle=1 si cycle négatif atteignable. */
Chemin bellman_ford_chemin(const Graphe* g, int src, int dst, int* has_negative_cycle);

/* ==========================
 * K plus courts chemins
 * ==========================
 * Implémentation de type Yen (simplifiée).
 * Retourne un tableau alloué de Chemin (à libérer par l’appelant via chemin_liberer sur chaque
 * entrée puis free du tableau). out_count = nombre réel de chemins (<=K).
 */
Chemin* k_plus_courts_chemins(const Graphe* g, int src, int dst, size_t K, size_t* out_count);

/* ==========================
 *   Redondance du réseau
 * ==========================
 * Indice [0,1] : 0 = disjoints, 1 = identiques.
 * consider_edges=1 : compare arêtes orientées ; sinon compare nœuds. */
float indice_redondance(const Chemin* chemins, size_t nb, int consider_edges);

#endif /* DIJKSTRA_H */
