# Brainstorming


## Bottleneck
### Cerca exhaustiva / arbre
Mirar totes les possibilitats, i per cada una executar els 2 algoritmes (baseline i el nostre)
+ Començar per una solució greedy

+ fer BFS per veure quins vertex estan més utilitzats pel baseline

+ calcular shortest paths i veure quins son més utilitzats pel baseline (tallar aquell putejar el baseline)

## Metaheurístiques
???



## Replanning

### Baseline (Shortest Path)
*It sorts the planning path requests of the same batch by weight from largest to smallest,
and then finds the global shortest path for each service in succession,
using only paths of invariant channels*


### Shortest path pel trencament
Es trenca el vertex _e_: _N_ → _M_

per tots els serveis que passen per _e_ (ordenats pel valor):
- busca el camí més curt de _N_ a _M_
- simplifica el camí si torna enrere pels mateixos vertexs per on ja passava el servei (si fa U-turns)

* Fer el shortest path fent un Dijkstra, amb costos a les arestes en funció de la capacitat restant

### Arbres
Cerca exhaustiva / arbres / branch & bound


