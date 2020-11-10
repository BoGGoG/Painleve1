Painleve I
========================================

Following the paper by Costin and Dunne [1904.11593](https://arxiv.org/abs/1904.11593).

ToDo
========================================
- Create the Series for the Painleve I
- Borel Transform


Creating the Series
----------------------------------------
$$y''(x) - 6 y^2(x) + x = 0$$

- want smooth real solution for large x
- Ecalle time $t = \frac{(24 x)^{5/4}}{30}$
- $y(x) = - \sqrt{\frac{x}{6}} ( 1 + h(t) )$
- now for $h(t)$ we have

$$\ddot{h}(t) + \frac{1}{t} \dot{h}(t) + h(t) \left( 1 + \frac{1}{2} h(t) \right)
- \frac{4}{25 t^2} \left(1 + h(t)\right) = 0$$
