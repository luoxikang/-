<center>原子物理课程报告</center>
=
<p align="right">罗熙康 湘雅医学院</p>
<p align="right">学号   2203150117</p>

一、简介
-
本课程报告选取的主题是“$\alpha$粒子的散射实验计算机模拟”我们选取的工具为Matlab.

这次报告里，我将会讲述自己的程序设计思路，程序的代码结构与解释，以及最后的结果及分析。

## 二、程序设计思路
### 1、$\alpha$粒子在单个金属原子电场中的运动轨迹
我们选取一个速度水平向，大小为$v$的$\alpha$粒子，从较远的地方开始飞向中心的Au原子。试着对其飞行轨迹进行建模

我们记录每个时刻的$\alpha$粒子的位置$\vec{x}_{t}$。Au原子核的位置我们假设其因为质量较大保持不动，记为$\vec{x}_{Au}$，因此我们有
$$
\vec{a}=\frac{Ze^2}{m_e} \cdot \frac{ (\vec{x}_{\alpha}-\vec x_{Au})}{|\vec{x}_{\alpha}-\vec x_{Au}|}\tag{1}
$$
其中
$$
\vec a=\frac{\mathrm{d}^2\vec{x}_\alpha}{\mathrm{d}t^2} 
$$

为了进行数值计算，我们对这个公式离散化$$
x_n=x_{n-1}+v_{x_n}\cdot dt \\
y_n=y_{n-1}+v_{y_n} \cdot dt
$$
其中$x_n,y_n$是在$t_n$时刻$\alpha$粒子的横纵坐标，$v_{x_n},v_{y_n}$是$\alpha$粒子在$t_n$时刻的速度的横纵方向上的分量  
此外还有

$$
v_{x_n}=v_{x_{n-1}}+a_{x_{n-1}} \cdot dt\\
v_{y_n}=v_{y_{n-1}}+a_{y_{n-1}} \cdot dt
$$
其中$a_{x_{n}},a_{y_n}$是在$t_n$时刻$\alpha$粒子加速度的横纵方向的分量。  
此外我们还有给出的初始值$x_0,y_0,a_{x_0},a_{y_0}$

按照以上公式进行不断的迭代计算，便可得出$\alpha$粒子的轨迹，并可进行后续的分析。

### 2. $\alpha$粒子在多个金属原子电场中的运动轨迹
我们假设多个Au原子核是按照下列方式分布的，其中原子核与核之间的间隔取金原子的直径$d=?pm$

我们只需要将公式$（1）$变为$$\vec{a_n}=\frac{Ze^2}{m_e} \cdot \frac{ (\vec{x}_{\alpha_n}-\vec x_{Au})}{|\vec{x}_{\alpha_n}-\vec x_{Au}|} \\ 
\vec{a}=\sum_{i=1}^{n}  \vec{a_n} \tag{1} $$
即可

## 三、程序的程序的代码结构与解释
`>源代码见附件`
整个程序的构架如下： 
### 基本参数
>L-金箔的长度
LBox-我们考虑的视野的大小
Z-金原子的电荷数
dt-我们数值计算时取的时间步长
RadiusAtom-金原子的半径
V0-金原子的初速度<
`注：为了计算及分析方便，我们参数取值与实际有不同`
```matlab
    L=10; %The Box  
    LBox=60;
    Z=10;%The Big Atom
    dt=0.001;%The Time Step
    RadiusAtom=10;%1.66e-8;%Atom's Diamter
    V0=10;

```
### 基本函数
* Initial_Alpha -初始化一个$\alpha$粒子的位置和初速度
* Initial_SingleBigAtom -在视野中央初始化单个Au原子
* Initial_MulBigAtom    -在视野中央初始化给定层数的Au原子
* Force -根据给定的$\alpha$粒子位置与Au原子位置计算$\alpha$粒子的加速度
* IfOutBox  -根据给定的$\alpha$粒子位置来判断其是否已经逃离视野
* Simulate  -根据给定的$\alpha$粒子的初始位置与初始速度，以及Au原子的位置，进行运动的模拟。其中我们假设Au因为自身质量原因是不动的。
* Gettheta  -根据给定的轨迹计算出速度的对应方向的角度
* ShowMotion    -根据给定的计算轨迹以及Au原子的位置画出运动示意图
* SingleAtom    -在中心有一个Au原子的情况下，根据给定的初始位置计算出最后的速度偏差角
* MulAtom   -大致同SingleAtom,区别在于此时需要给定Au原子的层数，对应的是多Au原子的状况。

```matlab
function Theta=SingleAtom(b,flag)
    Alpha=Initial_Alpha(b,V0);
    BigAtom=Initial_SingleBigAtom();
    AlphaMotion=Simulate(BigAtom,Alpha);
    if flag==1
    ShowMotion(AlphaMotion,BigAtom);
    end
    Theta=GetTheta(AlphaMotion);
end

function MulAtom(b,Layers)
    Alpha=Initial_Alpha(b,V0);
    BigAtom=Initial_MulBigAtom(Layers);
    AlphaMotion=Simulate(BigAtom,Alpha);
    ShowMotion(AlphaMotion,BigAtom);
end
    
function IfOut=IfOutBox(AlphaPos)
    x=AlphaPos(1); y=AlphaPos(2);
    if 0-LBox<=x && x<=LBox && 0-LBox<=y && y<=LBox
        IfOut=0;
    else
        IfOut=1;
    end
end

function AlphaMotion=Simulate(Ini_BigAtom,Initial_Alpha)
    IfOut=0;
    AlphaMotion=Initial_Alpha;
    while(~IfOut)
        AlphaPos=AlphaMotion(1:2,end); AlphaV=AlphaMotion(3:4,end);
        Force=Force(AlphaPos,Ini_BigAtom);
        AlphaV2=AlphaV+Force*dt;
        AlphaPos2=AlphaPos+AlphaV2*dt;
        AlphaMotion2=[AlphaPos2;AlphaV2];
        AlphaMotion=[AlphaMotion,AlphaMotion2];
        IfOut=IfOutBox(AlphaPos2);
    end     
end

function Ini_Alpha=Initial_Alpha(b,v)
    Alphax=-20; Alphay=0+b;
    Ini_Alpha=[Alphax;Alphay;v;0];
    
end

function Ini_BigAtom=Initial_SingleBigAtom()
    
    N0=5;% 5 Atoms each layer
    N=1;%total 
    Ini_BigAtom=zeros(2,N);
    NLay=floor(N/N0); NLastLay=N-N0*NLay; % The true number of layers is one more bigger than NLsy
    for i=1:NLay
        BigAtomx=0+(i-1)*RadiusAtom*2;
        for j=1:N0
            BigAtomy=(0+j)*L/(N0+1)-L/2;
            Ini_BigAtom(:,(i-1)*N0+j)=[BigAtomx;BigAtomy];
        end
    end
    BigAtomx=0+NLay*RadiusAtom*2;
    for j=1:NLastLay
        BigAtomy=(0+j)*L/(NLastLay+1)-L/2;
        Ini_BigAtom(:,N0*NLay+j)=[BigAtomx;BigAtomy];
    end
        
    
end


function Ini_BigAtom=Initial_MulBigAtom(Layers)
    
    N0=5;% 5 Atoms each layer
    N=Layers*N0;%total 
    Ini_BigAtom=zeros(2,N);
    NLay=floor(N/N0); NLastLay=N-N0*NLay; % The true number of layers is one more bigger than NLsy
    for i=1:NLay
        BigAtomx=0+(i-1)*RadiusAtom*2;
        for j=1:N0
            BigAtomy=(0+j)*L/(N0+1)-L/2;
            Ini_BigAtom(:,(i-1)*N0+j)=[BigAtomx;BigAtomy];
        end
    end
    BigAtomx=0+NLay*RadiusAtom*2;
    for j=1:NLastLay
        BigAtomy=(0+j)*L/(NLastLay+1)-L/2;
        Ini_BigAtom(:,N0*NLay+j)=[BigAtomx;BigAtomy];
    end
        
    
end

function Force=Force(Alpha,Ini_BigAtom)
    N=size(Ini_BigAtom,2);
    Force=[0;0]; AlphaPos=Alpha(1:2,:);
    for i=1:N
        BigAtomPos=Ini_BigAtom(:,i);
        RAlphaBigAtom=AlphaPos-BigAtomPos;
        Force0=Z/norm(RAlphaBigAtom)*RAlphaBigAtom;%*e^2/me;
        Force=Force+Force0;
    end
        
end

function Theta=GetTheta(AlphaMotion)
    AlphaVLast=AlphaMotion(3:4,end);
    TanTheta=AlphaVLast(2)/AlphaVLast(1);
    Theta=atan(TanTheta);
    Theta=cot(Theta/2);

end
    
function ShowMotion(AlphaMotion,Ini_BigAtom)
    figure();
    
    scatter(Ini_BigAtom(1,:),Ini_BigAtom(2,:),'filled','r');
    hold on
    plot(AlphaMotion(1,:),AlphaMotion(2,:));
    axis([-LBox,LBox,-LBox,LBox]);
end
```
`注：这些函数有和主函数共享变量，因此应该放置在主函数内`
### 主函数：单个Au原子下$\alpha$粒子的运动模拟
```matlab
  %% Single Atom
   b=10;
   SingleAtom(10,1)
```
结果如下

可以看到粒子的偏转角度十分之大

### 偏转角度与入射截距的关系
```matlab
 %% Show The Curve between b and Theta
    db=1;
    b=1:db:10;
    BTheta=zeros(2,length(b));
    for m=1:length(b)
        bi=b(m);
        Theta=SingleAtom(bi,0);
        BTheta(:,m)=[Theta;bi];
    end
    figure()
    plot(BTheta(1,:),BTheta(2,:));
```
结果如下
可以看出还是比较符合正性的线性关系的，图中X轴代表$\cot(\frac{\theta}{2})$,Y轴代表$b$
这与书上的公式
$$
b=\frac{a}{2} \cdot \cot \frac{\theta}{2}
$$
是相符的

###多Au原子下的运动
```matlab
%% Show The Motion In Mutiple BigAtom
    MulAtom(2,1);
```
可以看出来原子数量增多后，电场力也增大了，导致$\alpha$粒子反射角度更大了。

但我对此处存疑，可能是我的原子间距离选取的过小。还可以继续调节参数得出更详细的结果。

##三、总结与分析
+ 本次的模拟实验较成功，虽然多Au原子的情况参数还有调节的空间。
+ 偏转角度越大，意味着对应的横截距离约小
