
#ifndef BACKTRACKING_H
#define BACKTRACKING_H

#include "dijkstra.h"

/* =========================================
 *   Chemin avec contraintes (Backtracking)
 * =========================================
 * Objectif : minimiser la latence sous contraintes :
 *  - bande_passante_min  (>=)
 *  - cout_max            (<=)
 *  - must_visit[]        (passage obligatoire par ces nœuds)
 *  - exclude_nodes[]     (nœuds interdits)
 *  - appliquer_securite  (si 1 : n’utilise pas les arêtes non-sécurisées)
 *
 * Méthode : DFS backtracking avec élagage :
 *  - borne inférieure (LB) via Dijkstra sur le graphe renversé depuis dst
 *  - pruning si :
 *      * coût courant > cout_max
 *      * bande passante courante < bande_passante_min
 *      * latence courante + LB >= meilleure latence connue
 */

typedef struct Contraintes {
    const int* must_visit;
    size_t     must_visit_count;

    const int* exclude_nodes;
    size_t     exclude_count;

    float bande_passante_min;
    float cout_max;

    int   appliquer_securite;   /* 1 = interdire arêtes a->securite == 0 */
} Contraintes;

/* Calcule le meilleur chemin faisable (min latence). Retourne Chemin vide si impossible. */
Chemin chemin_contraint_backtracking(const Graphe* g, int src, int dst, const Contraintes* cs);

#endif /* BACKTRACKING_H */
