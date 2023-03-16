---
title: "Aproximación a optimización estructural"
date: 2023-02-16T11:48:52-05:00
draft: false
---

En este post se hablará de optimización estructural desde un punto más general, con el objetivo de afianzar las bases para futuras aplicaciones.

Para entender de forma simple el trasfondo de la optimización estructural es buena idea separar los dos términos y tomar la definición de cada uno. La optimización consiste en hacer las cosas mejor y una estructura se puede entender como un ensamble de material que tiene como propósito soportar cargar. Por lo tanto, la optimización estructural sé puede definir como la técnica de hacer que un ensamble de material soporte cargas de la mejor manera [1]. Esta definición es una idea abstracta de lo que es la optimización estructural, pero es un acercamiento simple que supone una idea más intuitiva del concepto. Cuando se menciona hacer lo “mejor” de algo por lo general lo que sé busca es minimizar y maximizar, un caso común en ingeniería es que se busque reducir los gastos de un proceso, esto se traduce a minimizar el costo de operación.

La optimización estructural utiliza un modelo matemático, este consiste en identificar una función objetivo, esta función es la que clasifica el diseño y será el su resultado el que se buscara minimizar o maximizar. Variable de diseño, es una función o vector que define el de sino, puede ser geometría o material. Variable de estado, representa la respuesta de la estructura, como los esfuerzos, deformaciones, desplazamientos, etc.

La optimización estructural se puede dividir en tres casos:

- Optimización de tamaño: es la minimización o maximizaron de las dimensiones del modelo.
- Optimización de forma: se refiere al dominio e.i en un problema de valores en la frontera, la optimización de forma se refiere a escoger el dominio que hará que la estructura presénteme menores esfuerzos.
- Optimización topológica: es la configuración de material, la optimización topológica permite la eliminación de material creando agujeros en el modelo, esto puede ser con el objetivo de minimizar la cantidad de material. Este caso es el que más se estudia en ingeniería.

Los desplazamientos, esfuerzos, deformación, peso, geometría, fuerzas, rigidez son algunas de las variables que pueden optimizarse. Un problema de optimización consiste en elegir alguna de estas variables como función a minimizar o maximizar y usar otra diferente como restricción.

### Referencias

1. Peter W. Christensen, A. K. (2009). An Introduction to Structural Optimization. Recuperado de https://link.springer.com/book/10.1007/978-1-4020-8666-3
