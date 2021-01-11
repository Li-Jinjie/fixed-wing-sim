## 记录mavsim值得我学习的地方

### Chapter 1 - 3

1. 所有可以重复使用的变量提前存在parameters中，比如各种gamma，在程序中只是简单的加和乘，减少重复计算过程。
2. 变量名取的和书上一样，取的真是太好了，完全可以像书上一样敲公式。
3. class内部的函数区分权限。_开头的函数和变量是半私有变量，在类外无法调用。这比我之前毕设时候的类写的好太多了。以后在写类的时候一定要区分好不同函数和参数的作用域。
4. 全用numpy可以提速，因为numpy底层是C实现的。在numpy中，用@表示矩阵乘法， 用[[],[]]的形式表示矩阵和列向量，将矩阵和列向量的形式统一起来。
5. 在第三章中，用四阶龙格库塔法求微分方程，真是活学活用！
6. 画的界面特别直观，还能拖动，用pyqt画的。

### Chapter 4