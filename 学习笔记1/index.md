# 学习记录(1)

## 感知机

![1701157495798](image/学习笔记(1)/1701157495798.png)

x是输入信号，y是输出信号，w(weight)是权重

输入信号分别乘以固定的权重(w*x)，传输至下一个神经元，超过某一个阈值(theta)时输出1，表示“神经元被激活”

![1701157752016](image/学习笔记(1)/1701157752016.png)


通过改变感知机的参数（w,theta）,也就是权重和阈值，相同构造的感知机可以表示与门、与非门、或门

![1701158276912](https://file+.vscode-resource.vscode-cdn.net/e%3A/my_website/content/posts/DLFS%281%29/image/%E5%AD%A6%E4%B9%A0%E7%AC%94%E8%AE%B0(1)/1701158276912.png "与非门真值表")![1701158333366](https://file+.vscode-resource.vscode-cdn.net/e%3A/my_website/content/posts/DLFS%281%29/image/%E5%AD%A6%E4%B9%A0%E7%AC%94%E8%AE%B0(1)/1701158333366.png "与门正值表")

![1701158300550](image/学习笔记(1)/1701158300550.png "或门真值表")


### 简单实现

```python
#与门
def AND(x1, x2):
    w1, w2, theta = 0.5, 0.5, 0.7
    tmp = w1*x1 + w2*x2
    if tmp <= theta:
        return 0
    elif tmp > theta:
        return 1
```

与非门的参数可以直接时 与门参数值的符号取反

### 导入权重与偏置

![1701158509935](image/学习笔记(1)/1701158509935.png)

将阈值（theta）换成-b,公式变形，也就是上式，同样可以表达感知机

感知机会计算输入信号和权重的乘积，然后加上偏置，通过与0大小判断来决定是否激活神经元

#### numpy实现感知机

```python
import numpy as np

def AND(x1, x2):
    x = np.array([0,1]) #input
    w = np.array([0.5, 0.5]) #weight
    b = -0.7 #bias
  
    tmp = np.sum(w*x) + b
    if tmp <= 0:
        return 0
    else:
        return 1

#直接修改w b可以分别实现与非门以及或门
```

w权重用于控制输入信号重要性占比，(-theta)偏置(b)用于调整神经元被激活的容易程度


### 局限性

单层感知机无法直接实现异或门（仅有一者为整输出1）（异或：拒绝其他）

![1701159459288](image/学习笔记(1)/1701159459288.png "感知机可视化")

上述构造的感知机，在坐标系中是通过一条直线划分成符合要求的两块空间，（以真值表中四个点为例）异或门只能用曲线分开

也就是单层感知机无法分离非线性空间


### 多层感知机

![1701159853667](image/学习笔记(1)/1701159853667.png)

![1701159865088](image/学习笔记(1)/1701159865088.png)

python实现异或门

```python
def AND(x1, x2):
    x = np.array([x1,x2]) #input
    w = np.array([0.5, 0.5]) #weight
    b = -0.7 #bias
  
    tmp = np.sum(w*x) + b
    if tmp <= 0:
        return 0
    else:
        return 1
  
def NAND(x1, x2):
    x = np.array([x1,x2]) #input
    w = np.array([-0.5, -0.5]) #weight
    b = 0.7 #bias
  
    tmp = np.sum(w*x) + b
    if tmp <= 0:
        return 0
    else:
        return 1
  
def OR(x1, x2):
    x = np.array([x1,x2]) #input
    w = np.array([0.5, 0.5]) #weight
    b = -0.2 #bias
  
    tmp = np.sum(w*x) + b
    if tmp <= 0:
        return 0
    else:
        return 1
  
def XOR(x1, x2):
    s1 = NAND(x1, x2)
    s2 = OR(x1, x2)
  
    return AND(s1, s2)

print ("0,0", XOR(0,0))
print ("1,0", XOR(1,0))
print ("0,1", XOR(0,1))
print ("1,1", XOR(1,1))
```

```python
0,0 0
1,0 1
0,1 1
1,1 0
```


### 与非门至计算机



## 神经网络

