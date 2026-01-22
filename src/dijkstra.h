
#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <stddef.h>
#include "graphe.h"  // Doit définir: Graphe, Arete, Noeud

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Représente un chemin calculé dans le graphe.
 */
typedef struct Chemin {
    int* noeuds;                 /**< Suite d'IDs de nœuds (source -> destination) */
    int  longueur;               /**< Nombre de nœuds dans le chemin */
    float latence_totale;        /**< Somme des latences le long du chemin */
    float cout_total;            /**< Somme des coûts le long du chemin */
    float bande_passante_min;    /**< Minimum de bande passante rencontrée (goulot) */
} Chemin;

/**
 * @brief Libère une structure Chemin allouée par les algorithmes.
 */
void chemin_liberer(Chemin* c);

/**
 * @brief Calcule le plus court chemin (par latence) avec Dijkstra.
 *
 * @param g           Graphe (orienté ou non) avec pondération Arete.latence >= 0
 * @param source      ID du nœud source (0..g->nb_noeuds-1)
 * @param destination ID du nœud destination
 * @param code_retour (optionnel) 0 si OK, 1 si source/destination invalides, 2 si pas de chemin,
 *                    3 si arête négative détectée (non supportée par Dijkstra)
 * @return Chemin*    Pointeur vers le chemin, ou NULL en cas d'échec (voir code_retour)
 *
 * @complexity O((V + E) * n_insertion) avec file de priorité sur liste chaînée.
 *             Ici, insertion triée O(n) et pop O(1) -> ~O(V^2 + E*logique) en pratique.
 */
Chemin* dijkstra_plus_court_par_latence(const Graphe* g, int source, int destination, int* code_retour);

/**
 * @brief Variante : remplit les tableaux des distances et parents, sans reconstruire le chemin.
 * @param dist   tableau de taille V (sortie), rempli avec la distance min par latence
 * @param parent tableau de taille V (sortie), parent[v] = prédécesseur de v dans l'APM, -1 pour source
 * @return 0 si OK, sinon != 0
 */
int dijkstra_remplir_tableaux(const Graphe* g, int source, float* dist, int* parent);

#ifdef __cplusplus
}
#endif

#endif /* DIJKSTRA_H */
