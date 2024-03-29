# Лабораторная работа № 6. Безусловный экстремум.

Выполнил студент группы №428  
Берсенев Даниил Романович 

## Вариант № 17
Задание: Найти точку **максимума** функции

![F(x1,x2,x3).jpg](Formuli/F(x1,x2,x3).jpg)

![[x1x2x3].jpg](Formuli/[x1x2x3].jpg)

методом **сопряженных градиентов**. Для одномерной минимизации использовать метод **квадратичной интерполяции**. Для поиска интервала унимодальности использовать алгоритм **скользящего окна**.
В окрестности точки максимума построить линии уровня и траекторию поиска (на одном графике).   


## Теоретическая часть

**Градиент** функции в точке-вектор, координаты которого-это частные производные по соответствующим аргументам(их значения в данной точке).

![grad.jpg](Formuli/grad.jpg)

Большинство процессов, используемых для приближенного решения задачи можно представить как итерационные в виде

![1.jpg](Formuli/1.jpg)

где ![2.jpg](Formuli/2.jpg) - вектор, который определяет направление движения от точки ![3.jpg](Formuli/3.jpg) к точке ![4.jpg](Formuli/4.jpg) - числовой множитель, величина которого задает длину шага в направлении ![2.jpg](Formuli/2.jpg).



**Метод сопряженних градиентов**

В этом методе ![ak.jpg](Formuli/ak.jpg) выбирается из условия 

![yclovie.jpg](Formuli/yclovie.jpg). 

При этом вектор ![2.jpg](Formuli/2.jpg) зависит не только от градиента функции ![fx^k.jpg](Formuli/fx^k.jpg), но и от градиента в предыдущей точке ![gradfx^k-1.jpg](Formuli/gradfx^k-1.jpg) (т.е. метод является двухшаговым), и строится либо по правилу 

**А)**![a1.jpg](Formuli/a1.jpg)

![a2.jpg](Formuli/a2.jpg)

либо по правилу 

**Б)**![b1.jpg](Formuli/b1.jpg)

![b2.jpg](Formuli/b2.jpg)

***n*** - размерность пространства независимых переменных. 
Вариант Б) отличается от варианта А) тем, что содержит так называемую процедуру обновления - для каждого ***k***, кратного ***n***, переход из точки ![3.jpg](Formuli/3.jpg) в точку ![x^k+1.jpg](Formuli/x^k+1.jpg) выполняется как в методе наискорейшего спуска.

**Метод квадратичной интерполяции**:

Здесь задаются три пробные точки ![x1=(a+b)na2.jpg](Formuli/x1=(a+b)na2.jpg), ![x2.jpg](Formuli/x2.jpg) и ![x3.jpg](Formuli/x3.jpg).

Для нахождения точки ![x2.jpg](Formuli/x2.jpg) задается шаг ![hbolshe0.jpg](Formuli/hbolshe0.jpg) в положительном направлении от точки ![x1.jpg](Formuli/x1.jpg), т.е. ![x2=x1+h.jpg](Formuli/x2=x1+h.jpg) и если ![f(x1)bolshef(x2).jpg](Formuli/f(x1)bolshef(x2).jpg) то ![x3=x1+2h.jpg](Formuli/x3=x1+2h.jpg), иначе ![x3=x1-h.jpg](Formuli/x3=x1-h.jpg).

Вычисляются значения функции в этих точках ![f(x1),f(x2),f(x3),.jpg](Formuli/f(x1),f(x2),f(x3),.jpg) строится квадратичный интерполяционый многочлен по трем точкам и находится его точка минимума по формуле

![x^.jpg](Formuli/x^.jpg).

Находится также точка ![xmin.jpg](Formuli/xmin.jpg).

Если знаменатель в формуле для нахождения минимума квадратичного интерполяционного многочлена равен нулю, т.е. все три точки лежат на одной прямой рекомендуется выбрать за ![x1=xmin.jpg](Formuli/x1=xmin.jpg) и повторить нахождение точки ![x^2.jpg](Formuli/x^2.jpg).

Критерием окончания этого процесса является выполнение условий для заданного ![epsilon.jpg](Formuli/epsilon.jpg)

![1.1.jpg](Formuli/1.1.jpg)

Если условия окончания не выполняются и

![x^3.jpg](Formuli/x^3.jpg)

точка ![x1.jpg](Formuli/x1.jpg) заменяется на точку ![argmin.jpg](Formuli/argmin.jpg), в противном случае точка ![x1.jpg](Formuli/x1.jpg) заменяется на ![x^2.jpg](Formuli/x^2.jpg).

**Алгоритм скользящего окна**

Для выбранной исходной точки ![x0.jpg](Formuli/x0.jpg) и выбраного окна шириной
![2hbolshe0.jpg](Formuli/2hbolshe0.jpg) около точки **x0** проверяeтся условие унимодальности

![yy.jpg](Formuli/yy.jpg)

Если условие выполнено то считается, что интервал унимодальности найден, в противном случае проверяется условие

![y.jpg](Formuli/y.jpg)

Если последнее выполнено, тогда окно сдвигается в право от точки ![x.jpg](Formuli/x.jpg) на ![hna2.jpg](Formuli/hna2.jpg), иными словами точка ![x.jpg](Formuli/x.jpg) изменяется на точку ![x=x+hna2.jpg](Formuli/x=x+hna2.jpg).

В противном случае окно сдвигается в лево от точки ![x.jpg](Formuli/x.jpg) на ![hna2.jpg](Formuli/hna2.jpg),
иными словами точка ![x.jpg](Formuli/x.jpg) изменяется на на точку ![x=x-hna2.jpg](Formuli/x=x-hna2.jpg).

Выбор ширины окна определяется экспериментально и целиком зависит от интуиции исследователя.

## Практическая часть
Программа содержит ряд ключевых *функций*:
* **Func** возвращает значение исходной функции в точке x.
* **Proizvx1** - производная по первой переменной.
* **Proizvx2** - производная по второй переменной.
* **Proizvx3** - производная по третьей переменной.
* **Tochka** - значение в точке *x+ap*.
* **SkolzOkoshko** - метод скользящего окна, возвращает середину отрезка унимодальности.
* **KvadrInt** - метод квадратичной интерполяции, возвращает точку максимума функции.
* **SoprGrad** - метод сопряженных градиентов, выводит точку максимума функции и значение функции в этой точке, а также создает файл ans1.dat и записывает в него точки траектории поиска.
Последовательность запуска: компиляция части программы на с++, ее запуск, затем запускаем часть на python, в которой находится отрисованная траектория поиска и линии уровня. 

### Результаты
В результате работы программы у функции 

![F(x1,x2,x3).jpg](Formuli/F(x1,x2,x3).jpg)

был найден экстремум в точке [6.8771,3.63493,-4.12432] (начальная точка [-2,-4,3]) за две итераций с точностью 1e-5. Ниже приведены рисуноки с изображением линий уровня анализируемой функции и траектория поиска экстремума:
![x1x2.png](Formuli/x1x2.png)
![x1x3.png](Formuli/x1x3.png)
![x2x3.png](Formuli/x2x3.png)
