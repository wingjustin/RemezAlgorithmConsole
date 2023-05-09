# RemezAlgorithmConsole
Remez Algorithm

<h3>雷米兹交换演算法</h3>

<p>
用切比雪夫插值作為初始插值多項式
<br><br>
切比雪夫插值可提高對插值誤差的控制, n+1個插值點的多項式插值誤差的分子為一元n+1階切比雪夫多項式乘函數f(c)的n+1階導數,分母為n+1階乘,x在[-1,1], c在最少插值點和最大插值點之間
<br><br>
一元n+1階切比雪夫多項式的絕對最大值為1/(2^n), 所以多項式插值絕對誤差最大值控制在f^(n+1)(c)/(2^n)/(n+1)!內
<br><br>
而且對其它任意一元n+1階函數g(x),x在[-1,1], 其絕對最大值必定大於等於1/(2^n)
<br><br>
</p>

<div>References</div>
<ul>
<li>1. https://youtu.be/vGTZVerkTds?t=2212</li>
<li>2. https://zhuanlan.zhihu.com/p/393903459</li>
<li>3. https://blog.csdn.net/weixin_42664622/article/details/103673875</li>
</ul>
<hr>
<p>
然後找最大或最小絕對誤差或相對誤差值作為控制點,再解方程,不斷循環
<br><br>
,過程主要用到TOMS748 求根算法和Brent's Method 求極值算法, 兩種算法都不需要求對函數f(x)求導, 加上LU分解法求解
<br><br>
  支持有理函數逼近及權重函數。當控制點正負沒有交錯會自動嘗試找尋附近相反符號點作為新控制點。
  <br><br>
</p>
<div>Parameters</div>
<ul>
<br>
<li>func = 目標函數f(x)</li>
<br>
<li>weight_func = 權重函數w(x)</li>
<br>
<li>oN = 分子階數</li>
<br>
<li>oD = 分母階數,階數為零時為逼近函數為多項式</li>
<br>
<li>pinned [true|false] = true 為 控制逼近函數必須通過原點(0,0), 並且盡量避免控制點為0導致出現除零情況</li>
<br>
<li>skew [-100,100] 為使初始切比雪夫零點左偏或右偏, 正數為左偏, 負數為右偏, 預設為0</li>
<br>
<li>relative_error[true|false] = 如果沒有權重函數, true 為 對相對誤差最佳化, false 為 對絕對誤差最佳化, 預設為true。 當有權重函數時此設定無意義。</li>

</ul>

<hr>

<div>References</div>
<ul>
<li>1. https://www.boost.org/doc/libs/1_46_1/libs/math/doc/sf_and_dist/html/math_toolkit/backgrounders/remez.html</li>
<li>2. https://www.boost.org/doc/libs/1_82_0/libs/math/doc/html/math_toolkit/internals/minimax.html</li>
<li>3. boost.org, .\libs\math\include_private\boost\math\tools\remez.hpp</li>
<li>4. https://github.com/samhocevar/lolremez/wiki</li>
<li>5. https://github.com/samhocevar/lolremez/blob/main/src/solver.cpp</li>
</ul>
