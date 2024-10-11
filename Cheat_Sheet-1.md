# Cheat Sheet: Optimización de Continuidad del Servicio en Redes Ópticas

---

## Definiciones Clave

- **Red Óptica**: Se modela como un grafo no dirigido y conectado. Los **nodos** representan routers y las **arestas** representan fibras ópticas.
- **Servicio Óptico**: Es una conexión entre un nodo origen y uno destino, representado como un camino en el grafo. Cada servicio tiene un **ancho de banda** (cantidad de longitudes de onda necesarias), un **valor de servicio** y un **camino**.
- **Longitudes de Onda**: Son las frecuencias disponibles en cada fibra (aresta). Los servicios usan canales, que consisten en un conjunto contiguo de longitudes de onda.
- **Conversión de Canales**: Oportunidad de cambiar la longitud de onda entre dos arestas en un nodo. Los nodos tienen un número limitado de estas oportunidades.
- **Fallo en la Red (Fiber Cut)**: Ocurre cuando una fibra óptica (aresta) se daña, lo que obliga a los servicios que la usaban a ser replanteados.

---

## Nomenclatura

- **N**: Número de nodos en la red.
- **M**: Número de arestas (fibras ópticas) en la red.
- **S**: Número de servicios iniciales que se ejecutan en la red.
- **Wavelength**: Longitud de onda, se refiere a la frecuencia utilizada en una aresta.
- **Bandwidth (B)**: Número de longitudes de onda contiguas necesarias para un servicio.
- **Service Path**: Camino que sigue un servicio desde el nodo origen al nodo destino.
- **Channel**: Un conjunto de longitudes de onda que son utilizadas por un servicio en cada aresta de su camino.
- **Conversion Opportunity**: Número de oportunidades que tiene un nodo para cambiar la longitud de onda entre dos arestas.
- **Service Value (V)**: Valor de un servicio, se busca maximizar el valor total de los servicios que sobreviven después de varios fallos.
- **Edge Failure**: Indica un fallo en una fibra óptica (aresta), lo que requiere replanificar los servicios que la utilizaban.

---

## Restricciones

1. **Servicios sin Repetición de Nodos o Arestas**:
   - Un servicio óptico se representa como un camino **simple** (sin nodos ni arestas repetidos) desde el nodo origen al nodo destino.

2. **Uso de Longitudes de Onda**:
   - Para cada servicio, las longitudes de onda utilizadas en cada aresta del camino deben ser **contiguas**.
   - Diferentes servicios pueden compartir la misma aresta, pero no pueden usar la misma longitud de onda en esa aresta.

3. **Conversión de Canales**:
   - Los nodos pueden cambiar la longitud de onda entre dos arestas mediante las oportunidades de **conversión de canales**.
   - Cada nodo tiene un número limitado de conversiones disponibles. Una vez que se agotan, no se pueden volver a utilizar.
   - Un servicio debe usar el **mismo conjunto de longitudes de onda en todas las arestas** de su camino, a menos que se utilice una conversión de canal en un nodo intermedio. **Cada vez que se hace una conversión de canal, cuenta como una oportunidad consumida**.

4. **Fallos en la Red**:
   - Si una fibra óptica (aresta) falla, los servicios que usan esa aresta deben ser **replanificados**.
   - La replanificación debe cumplir las mismas restricciones de origen, destino y ancho de banda que el servicio original.
   - Un servicio puede reutilizar sus propios recursos previos (longitudes de onda y conversiones), pero **no puede usar recursos de otros servicios**.

5. **Servicios Muertos**:
   - Si no es posible replantear un servicio, se considera **muerto**. Los recursos utilizados por ese servicio en su camino original no se liberarán y no podrá ser replanificado más adelante.

6. **Planificación Simultánea**:
   - Todos los servicios afectados por un fallo deben ser replanteados simultáneamente. Los recursos del camino original de los servicios solo se liberan una vez que la replanificación ha finalizado.

---

## Interacción y Salida Esperada

1. **Secuencias de Fallos**:
   - Debes proporcionar **grupos de secuencias de fallos** de fibras ópticas (arestas) para probar la red. Cada secuencia de fallos debe ser única y no tener alta similitud entre ellas (medido por el **Índice de Jaccard**).
   - La secuencia de fallos se especifica por el número de arestas que fallan y sus identificadores.

2. **Replanificación de Servicios**:
   - Después de recibir un fallo, debes replantear los servicios afectados y proporcionar el número de servicios replanificados con éxito.
   - Para cada servicio replanificado, debes especificar el nuevo camino (secuencia de arestas) y las longitudes de onda utilizadas en cada aresta.

---

## Criterios de Puntaje

- **Valor de Servicio**: El puntaje se basa en el valor total de los servicios que logran sobrevivir tras los fallos.
- **Comparación con Solución Baseline**: Tu solución se compara con una **solución baseline** que utiliza el algoritmo de camino más corto y selecciona caminos en función del peso de los servicios.
- **Maximización**: Debes maximizar el puntaje total tratando de replantear la mayor cantidad de servicios de valor alto.

---

## Errores Comunes

1. **Número de Servicios Incorrecto**: Error en la cantidad de servicios replanteados.
2. **Duplicación de IDs**: Servicios o arestas repetidos en la salida.
3. **Anchura de Servicio Inconsistente**: Uso incorrecto del ancho de banda (longitudes de onda) en las arestas.
4. **Camino Desconectado**: El camino replanteado no conecta correctamente los nodos origen y destino.
5. **Insuficiencia de Oportunidades de Conversión**: Se agotan las oportunidades de conversión de canal en los nodos.

---

## ¿Por qué Proporcionar las Arestas Fallidas?

En este desafío, se te pide generar **tus propias secuencias de fallos**. Esto puede estar relacionado con la evaluación del algoritmo en diferentes escenarios. Al dar tus propias secuencias de fallos, el objetivo es que tu algoritmo no solo sea bueno para los fallos predefinidos por el sistema, sino que también logres identificar casos problemáticos o "cuellos de botella" que afecten al rendimiento de la solución baseline. Así, se prueba la robustez de tu algoritmo ante diferentes configuraciones de fallos.

Tu tarea es doble:
1. **Generar secuencias de fallos** que puedan poner a prueba la red.
2. **Replanificar** los servicios para maximizar la cantidad de servicios que sobreviven a estos fallos.

De esta forma, el sistema evalúa cómo maneja tu solución las secuencias que tú mismo propones.
