
#include "backtracking.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define INFINI_FLT (1e30f)

typedef struct {
    int*   pile;              /* nœuds courants */
    int    profondeur;        /* taille courante de la pile */
    int    capacite;          /* capacité alloc */
    float  latence_courante;
    float  cout_courant;
    float  bpmin_courante;    /* min(bande passante) le long du chemin courant */
} EtatParcours;

typedef struct {
    float meilleur_latence;
    int*  meilleur_chemin;
    int   meilleur_len;
    float meilleur_cout;
    float meilleur_bpmin;
} Meilleur;

/* =========================
 *  Helpers
 * ========================= */

static int contient(const int* tab, int n, int x) {
    if (!tab) return 0;
    for (int i = 0; i < n; ++i) if (tab[i] == x) return 1;
    return 0;
}

/* renvoie l'index de x dans 'obligatoires' ou -1 */
static int index_obligatoire(const int* obligatoires, int n, int x) {
    if (!obligatoires) return -1;
    for (int i = 0; i < n; ++i) if (obligatoires[i] == x) return i;
    return -1;
}

static const Arete* trouver_arete(const Graphe* g, int u, int v) {
    const Arete* a = g->noeuds[u].aretes;
    while (a) { if (a->destination == v) return a; a = a->suivant; }
    return NULL;
}

static void meilleur_init(Meilleur* m) {
    m->meilleur_latence = INFINI_FLT;
    m->meilleur_chemin = NULL;
    m->meilleur_len = 0;
    m->meilleur_cout = 0.0f;
    m->meilleur_bpmin = 0.0f;
}

static void meilleur_maj(Meilleur* m, const EtatParcours* etat) {
    free(m->meilleur_chemin);
    m->meilleur_len = etat->profondeur;
    m->meilleur_chemin = (int*)malloc(sizeof(int) * m->meilleur_len);
    if (m->meilleur_chemin) {
        memcpy(m->meilleur_chemin, etat->pile, sizeof(int) * m->meilleur_len);
    }
    m->meilleur_latence = etat->latence_courante;
    m->meilleur_cout    = etat->cout_courant;
    m->meilleur_bpmin   = (etat->bpmin_courante == INFINI_FLT ? 0.0f : etat->bpmin_courante);
}

static void etat_push(EtatParcours* etat, int v) {
    if (etat->profondeur >= etat->capacite) return; /* protection */
    etat->pile[etat->profondeur++] = v;
}

static void etat_pop(EtatParcours* etat) {
    if (etat->profondeur > 0) etat->profondeur--;
}

/* =========================
 *  Backtracking
 * ========================= */

static void dfs_contraint(const Graphe* g,
                          int current, int destination,
                          const ContraintesChemin* cons,
                          int* visites,
                          int* oblig_visites_flags, int* nb_oblig_couverts,
                          EtatParcours* etat,
                          Meilleur* meilleur) {
    /* Élagage par borne: si déjà pire que le meilleur, on coupe */
    if (etat->latence_courante >= meilleur->meilleur_latence) {
        return;
    }

    if (current == destination) {
        /* Vérifier que tous les passages obligatoires ont été rencontrés */
        if (cons->nb_passages_obligatoires > 0) {
            if (*nb_oblig_couverts < cons->nb_passages_obligatoires) {
                return; /* pas tous couverts */
            }
        }
        /* Chemin valide et optimal mieux que meilleur courant */
        meilleur_maj(meilleur, etat);
        return;
    }

    /* Explorer les successeurs */
    for (const Arete* e = g->noeuds[current].aretes; e; e = e->suivant) {
        int v = e->destination;
        if (v < 0 || v >= g->nb_noeuds) continue;

        /* Nœud exclu ? */
        if (contient(cons->exclus, cons->nb_exclus, v)) continue;

        /* Éviter cycles */
        if (visites[v]) continue;

        /* Contraintes d'arête */
        float new_bpmin = etat->bpmin_courante;
        if (new_bpmin == INFINI_FLT) new_bpmin = e->bande_passante;
        else if (e->bande_passante < new_bpmin) new_bpmin = e->bande_passante;

        if (cons->bande_passante_min > 0.0f && new_bpmin < cons->bande_passante_min) {
            continue; /* arête trop étroite */
        }

        float new_cout    = etat->cout_courant + e->cout;
        if (new_cout > cons->cout_max) {
            continue; /* dépassement coût max */
        }

        float new_latence = etat->latence_courante + e->latence;

        /* Élagage par meilleure latence connue */
        if (new_latence >= meilleur->meilleur_latence) {
            continue;
        }

        /* Marquages */
        visites[v] = 1;
        int idx_ob = index_obligatoire(cons->passages_obligatoires, cons->nb_passages_obligatoires, v);
        int added_oblig = 0;
        if (idx_ob >= 0 && oblig_visites_flags && !oblig_visites_flags[idx_ob]) {
            oblig_visites_flags[idx_ob] = 1;
            (*nb_oblig_couverts)++;
            added_oblig = 1;
        }

        /* Empiler et maj métriques */
        etat_push(etat, v);
        float old_lat  = etat->latence_courante;
        float old_cout = etat->cout_courant;
        float old_bp   = etat->bpmin_courante;

        etat->latence_courante   = new_latence;
        etat->cout_courant       = new_cout;
        etat->bpmin_courante     = new_bpmin;

        dfs_contraint(g, v, destination, cons, visites, oblig_visites_flags, nb_oblig_couverts, etat, meilleur);

        /* Désempiler et restauration */
        etat->latence_courante = old_lat;
        etat->cout_courant     = old_cout;
        etat->bpmin_courante   = old_bp;
        etat_pop(etat);

        if (added_oblig) {
            oblig_visites_flags[idx_ob] = 0;
            (*nb_oblig_couverts)--;
        }
        visites[v] = 0;
    }
}

Chemin* backtracking_chemin_contraint(const Graphe* g, int source, int destination,
                                      const ContraintesChemin* cons, int* code_retour) {
    if (code_retour) *code_retour = 0;
    if (!g || !cons) { if (code_retour) *code_retour = 1; return NULL; }
    if (source < 0 || destination < 0 || source >= g->nb_noeuds || destination >= g->nb_noeuds) {
        if (code_retour) *code_retour = 1;
        return NULL;
    }
    if (contient(cons->exclus, cons->nb_exclus, source) ||
        contient(cons->exclus, cons->nb_exclus, destination)) {
        if (code_retour) *code_retour = 2;
        return NULL;
    }

    /* État initial */
    int V = g->nb_noeuds;
    int* visites = (int*)calloc((size_t)V, sizeof(int));
    if (!visites) { if (code_retour) *code_retour = 1; return NULL; }

    int* oblig_flags = NULL;
    int  nb_oblig_couverts = 0;
    if (cons->nb_passages_obligatoires > 0) {
        oblig_flags = (int*)calloc((size_t)cons->nb_passages_obligatoires, sizeof(int));
        if (!oblig_flags) { free(visites); if (code_retour) *code_retour = 1; return NULL; }
        /* Si la source est déjà obligatoire, marquer */
        int idx_src = index_obligatoire(cons->passages_obligatoires, cons->nb_passages_obligatoires, source);
        if (idx_src >= 0) { oblig_flags[idx_src] = 1; nb_oblig_couverts = 1; }
    }

    EtatParcours etat = {0};
    etat.capacite = V;
    etat.pile = (int*)malloc(sizeof(int) * (size_t)V);
    if (!etat.pile) {
        free(visites); free(oblig_flags);
        if (code_retour) *code_retour = 1;
        return NULL;
    }
    etat.profondeur = 0;
    etat.latence_courante = 0.0f;
    etat.cout_courant = 0.0f;
    etat.bpmin_courante = INFINI_FLT;

    Meilleur meilleur; meilleur_init(&meilleur);

    /* Option d'élagage fort : on peut initialiser meilleur.meilleur_latence
       avec Dijkstra (latence pure), puis chercher un meilleur sous contraintes. */
    int rc_dij = 0;
    Chemin* c0 = dijkstra_plus_court_par_latence(g, source, destination, &rc_dij);
    if (c0) {
        meilleur.meilleur_latence = c0->latence_totale; /* borne supérieure initiale */
        chemin_liberer(c0);
    }

    /* Départ */
    visites[source] = 1;
    etat_push(&etat, source);
    int idx_src_ob = index_obligatoire(cons->passages_obligatoires, cons->nb_passages_obligatoires, source);
    if (idx_src_ob >= 0 && oblig_flags && !oblig_flags[idx_src_ob]) {
        oblig_flags[idx_src_ob] = 1;
        nb_oblig_couverts++;
    }

    dfs_contraint(g, source, destination, cons, visites, oblig_flags, &nb_oblig_couverts, &etat, &meilleur);

    /* Construction résultat */
    Chemin* res = NULL;
    if (meilleur.meilleur_chemin && meilleur.meilleur_len > 0 && meilleur.meilleur_latence < INFINI_FLT/2) {
        res = (Chemin*)malloc(sizeof(Chemin));
        if (res) {
            res->noeuds = meilleur.meilleur_chemin; /* transfert de propriété */
            res->longueur = meilleur.meilleur_len;
            res->latence_totale = meilleur.meilleur_latence;
            res->cout_total = meilleur.meilleur_cout;
            res->bande_passante_min = meilleur.meilleur_bpmin;
            /* Ne pas free meilleur.meilleur_chemin maintenant */
            meilleur.meilleur_chemin = NULL;
        } else {
            free(meilleur.meilleur_chemin);
        }
    } else {
        if (code_retour) *code_retour = 2; /* aucune solution */
        free(meilleur.meilleur_chemin);
    }

    /* Nettoyage */
    free(visites);
    free(oblig_flags);
    free(etat.pile);

    return res;
}
``
