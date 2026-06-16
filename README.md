Pregunta 1: Eficiencia Computacional: Tomando en cuenta los resultados de la terminal, explique la relación entre el número de evaluaciones de la función f(x) y el
error absoluto obtenido por Gauss-Legendre frente a los métodos de Newton-Cotes (Trapecio y Simpson). ¿Por qué se dice que Gauss optimiza el costo computacional?

Para que Trapecio y Simpson bajen su error, tienen que partir el área en un montón de pedazos (por eso hacen 13 evaluaciones). Gauss optimiza el costo computacional
porque no necesita tantos puntos; al hacer que la posición de los nodos no sea fija, se ajusta mucho mejor a la forma de la función, evalúa solo 5 veces,
ahorrando procesamiento, y encima da un mejor resultado.

Pregunta 2: El Efecto de la Oscilación: Observe el comportamiento de la función en el intervalo dado. ¿Qué ocurre con la pendiente y la oscilación en las zonas 
de mayor curvatura? Explique matemáticamente por qué los métodos compuestos tradicionales fallan más en esa zona que la cuadratura de Gauss.

En las partes donde la curva oscila mucho o cambia de dirección rápido (como en el coseno), los métodos tradicionales fallan porque sus puntos de evaluación están 
fijos, si la distancia entre puntos es muy grande, se terminan saltando los picos o valles de la gráfica, la cuadratura de Gauss no usa distancias iguales, 
sino que concentra sus puntos de forma estratégica (más hacia los bordes), lo que le permite no perderse esos cambios bruscos de la pendiente.

Pregunta 3: Límites del Método: Si cambiáramos la función por el polinomio exacto P(x) = 7x^7 - 3x^4 + 2x, ¿cuántos nodos n requeriría Gauss-Legendre para obtener
un error absoluto de exactamente cero? Sugerencia: Justifique utilizando el teorema del grado de precisión máxima.

Para que Gauss-Legendre dé un error de exactamente cero con el polinomio P(x) = 7x^7 - 3x^4 + 2x, usamos el teorema de su grado máximo de precisión: 
2n - 1 \ge d
Como nuestro polinomio es de grado 7, el grado es $d = 7$. Reemplazando en la inecuación:
2n - 1 \ge 7
2n \ge 8
n \ge 4
Se necesitan como minimo 4 nodos para obtener la solución exacta.
