
#ifndef BACKTRACKING_H
#define BACKTRACKING_H

#include <stdbool.h>
#include "graphe.h"
#include "dijkstra.h"  /* pour le type Chemin et chemin_liberer() */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Contraintes pour la recherche de chemin.
 */
typedef struct ContraintesChemin {
    float bande_passante_min;     /**< Bande passante minimale exigée sur chaque arête (Mbps). 0 si aucune. */
    float cout_max;               /**< Coût total maximal autorisé (<=). FLT_MAX si aucune limite. */

    const int* passages_obligatoires; /**< Tableau d'IDs de nœuds par lesquels le chemin doit passer (peut être NULL) */
    int        nb_passages_obligatoires;

    const int* exclus;            /**< Tableau d'IDs de nœuds interdits (peut être NULL) */
    int        nb_exclus;
} ContraintesChemin;

/**
 * @brief Recherche d'un chemin minimisant la latence sous contraintes multiples (backtracking + élagage).
 *
 * @param g           Graphe
 * @param source      ID source
 * @param destination ID destination
 * @param cons        Contraintes (voir structure)
 * @param code_retour (optionnel) 0=OK, 1=paramètres invalides, 2=aucune solution
 * @return Chemin*    Chemin optimal trouvé en latence, NULL si échec (voir code_retour)
 *
 * @complexity Exponentielle dans le pire cas O(b^d), avec élagages significatifs.
 */
Chemin* backtracking_chemin_contraint(const Graphe* g, int source, int destination,
                                      const ContraintesChemin* cons, int* code_retour);

#ifdef __cplusplus
}
#endif

#endif /* BACKTRACKING_H */
``
