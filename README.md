# Résilience des Essaims de Nano-Satellites

Projet académique (ENSEEIHT – 2A SN Télécoms) portant sur la modélisation et l'analyse de la résilience d'un essaim de 100 nano-satellites en orbite lunaire, dans le contexte d'un interféromètre radio spatial distribué (inspiré du projet NOIRE).

## Objectifs du projet
L'objectif est de modéliser mathématiquement la dynamique de l'essaim à l'aide de la théorie des graphes temporels afin d'évaluer sa robustesse et sa résilience face à des scénarios de pannes (épuisement énergétique, pannes franches).
- Modéliser l'essaim comme un graphe temporel non orienté et pondéré `Gt = (Vt, Et)` à partir de traces de mobilité orbitales.
- Évaluer la robustesse et la résilience du réseau via des métriques issues de la théorie des graphes.
- Simuler des scénarios de pannes aléatoires et ciblées pour mesurer l'impact sur la survie du système.

## Fonctionnalités
- **Construction du graphe** : Traitement des coordonnées (x, y, z) avec une portée adaptative à 3 paliers (20 km, 40 km, 60 km), traduite en poids de routage.
- **Métriques de robustesse** : Robustesse de flux, coût de routage global, efficacité réseau globale.
- **Métriques de résilience** : Redondance des chemins, disparité (différence symétrique des nœuds internes), criticité (betweenness centrality avec seuil).
- **Scénarios de suppression** : Simulation de pannes aléatoires, basées sur la centralité de degré, ou ciblant l'importance des nœuds (relais majeurs).
- **Visualisation** : Génération de topologies 3D interactives (Plotly) et de courbes d'évolution temporelle des métriques (Matplotlib).

## Résultats clés
- Le réseau nominal est quasi-connexe en permanence (robustesse de flux proche de 100 %), mais présente un coût de routage très élevé et une forte criticité structurelle (présence de goulets d'étranglement).
- Les suppressions ciblées ne fragmentent pas immédiatement le réseau (grâce à sa densité), mais détruisent les routes de secours (chute de la redondance et de la disparité).
- L'analyse met en évidence un effet de cascade : la perte des nœuds centraux transfère la charge de trafic aux nœuds survivants, faisant augmenter la criticité globale.
- Conclusion : L'exploitation de la topologie sous forme de graphe non divisé mène à un effondrement structurel. L'intégration de stratégies de division équitable du graphe (Fair Graph Division) est nécessaire.

## Stack technique
- **Python 3**
- **Bibliothèques** : `numpy`, `pandas`, `networkx`, `matplotlib`, `plotly`, `tqdm`.

## Structure du dépôt
- `simu.py` : Pipeline principal (importation des données, calcul des métriques globales, scénarios de pannes, génération des graphiques).
- `swarm_sim.py` : Modèle objet (classes `Node` et `Swarm`, calcul des métriques locales, algorithmes d'échantillonnage).
- `rapport_science_de_reseaux_version_finale.pdf` : Rapport complet du projet.
- `index.html` : Page web du projet.

## Installation et exécution
```bash
pip install numpy pandas matplotlib networkx plotly tqdm
python simu.py
