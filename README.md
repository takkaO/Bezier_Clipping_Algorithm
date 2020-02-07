# Bezier Clipping Algorithm

## Features
Implemented Bezier Clipping algorithm with Python.  
This program finds the intersection of a straight line and a Bezier curve. 

## Environment
- Python3.6
  - scipy == 1.3.1
  - numpy == 1.17.3
  - matplotlib == 3.1.1

## Run sample
Execute below command.  
```
python main.py
```

You will get below images. Red point is intersection.
![result1](https://github.com/takkaO/Bezier_Clipping_Algorithm/blob/images/Figure_1.png?raw=true)
```
Result1:
t:0.074 -> (x, y) = [ 0.09929439 -1.76824575]
t:0.595 -> (x, y) = [0.88230049 0.05870114]
t:0.878 -> (x, y) = [1.33648382 1.11846224]
```

![result2](https://github.com/takkaO/Bezier_Clipping_Algorithm/blob/images/Figure_2.png?raw=true)
```
t:0.067 -> (x, y) = [ 3.27081548 -3.36609438]
t:0.829 -> (x, y) = [7.00386554 1.4335414 ]
```


## Reference
- [Bezier Clipping](http://nishitalab.org/user/nis/ourworks/BezClip/BezierClipping.html)
- [Bezier Clipping およびその計算機支援形状設計への応用](https://www.ieice.org/jpn/event/FIT/2016/data/pdf/I-011.pdf)
- [Curve intersection using Bezier clipping](http://nishitalab.org/user/nis/cdrom/cad/CAGD90Curve.pdf)
- [曲線・曲面の処理](https://www.jstage.jst.go.jp/article/bjsiam/13/3/13_KJ00003509912/_pdf)
- [BEZIER CURVES INTERSECTION USING MONGE PROJECTION](http://www.sccg.sk/~kg/palaj/pub/3_Kocovce%202007.pdf)
- [ベジェ曲線を手で描いてみる - にせねこメモ](https://nixeneko.hatenablog.com/entry/2015/06/26/075022)
- [A Primer on Bezier Curves](https://pomax.github.io/bezierinfo/#splitting)
- [A Matrix method for efficient computation of Bernstein coefficients](https://interval.louisiana.edu/reliable-computing-journal/volume-17/reliable-computing-17-pp-40-71.pdf)