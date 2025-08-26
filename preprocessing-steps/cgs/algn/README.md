# algn

sqrt(2) を代数的数として扱った計算を記録

## Computing instructons

### Computing CGS

1. Programs in this directory use CGS program by Prof. Katsusuke Nabeshima.
1. Computation of the CGS of the basis F w.r.t. Lex order can be executed as follows:
    ```
    % asir
    load("cgs-F-input-algn-2.rr")$
    ```
    The output is stored in ```G2```.
    If the variable ```OUTPUT = 1```, then the contents of ```G2``` is stored in ```G2.dat```.
1. Or, w.r.t. DegRevLex order, do the following:
    ```
    % asir
    load("cgs-F-input-algn-0.rr")$
    ```
    The output is stored in ```G0```.
    If the variable ```OUTPUT = 1```, then the contents of ```G0``` is stored in ```G0.dat```.
1. Reverse the order of the elements in 
    ```
    % asir
    G2 = bload("G2.dat")$
    G2r = reverse(G2)$
    bsave(G2r, "G2-reverse.rr")$
    ```
### Extracting CGS

Run the script as follows:
````
% asir
load("extract-cgs-G-script.rr")$
````
Then, it reads the contents of ```G2-reverse.dat```. 
For the segment (S_i,G_i) with S_i = V_R(I_{i,1})\ V_R(I_{i,2}), its elements are extracted as <br />
* I_{i,1}: ```F-segments/F-i-1.rr```
* I_{i,2}: ```F-segments/F-i-2.rr```
* G_i: ```G-basis/G-i.rr```

### Verifying the segments

The presence of real points in each segment was verified by manual calculation.
The results of the verification for each segment S_i are recorded as follows:
* If S_i has real point(s): ```F-segments/F-i-verification.log``` 
* Otherwise:  ```F-segments/F-i-verification-false.log```

### Arrangements on the elements in G (Algorithm 2)

1. Arrange the elements in G as shown in Algorithm 2 as:
    ```
    % asir
    load("generate-G2-new.rr")$
    ```
    Then, the program reads ```G2-reverse.dat``` and output to the variable ```G2new```. if the variable ```OUTPUT = 1```, the contents of ```G2new``` is saved in ```G2-new.dat```. Note: ```arrange-G.rr``` substitutes ```a``` with ```2^(1/2)```.
1. ```generate-G2-new-a.rr``` can also be used for the same purpose. 
The difference between ```arrange-G.rr``` is that, with the use of ```arrange-G-a.rr```, the variable of ```a``` is preserved.

### Verifying zero-dimensional ideals

1. Verify zero-dimensional ideals as:
    ```
    % asir
    load("zero-dimensional-test-script-G2-new.rr")$
    ```
    Then, variables which do not appear in the leading monomial in the Groebner basis appears.

1. By using a modified program, you can check the leading coefficient in each element in the Groebner basis.
    ```
    % asir
    load("zero-dimensional-test-script-2-G2-new.rr")$
    ```

### Re-arranging the CGS

According to the above observation, re-arrange CGS as:

```
% asir
load("generate-G2-2.rr$)$
```
With this operation, the program reads ```G2-new.dat``` and ```G2-new-a.dat``` and outputs to the variables ```G2_2``` and ```G2_2a```, respectively. 
In this program, we reduce the number of CGS bases and adjust them so that each basis consists of six polynomials whose leading variables are arranged in the order [s_7, c_7, s_4, c_4, s_1, c_1].
If the variable ```OUTPUT = 1```, the contents are saved in ```G2-2.dat``` and ```G2-2-a.dat```, respectively.
Note: while in ```G2-2.dat```, the variable ```a``` is substituted with ```2^(1/2)```, in ```G2-2-a.dat```, the variable ```a``` remains.

For the subsequent calculation steps, refer to [[../../hermite/algn/README.md]].


## スクリプト

* cgs-F-input-algn.0.rr: 単項式順序 DegRevLex でCGSを計算するスクリプト。出力先はG0.dat。
* cgs-F-input-algn.2.rr: 単項式順序 Lex でCGSを計算するスクリプト。出力先はG2.dat。
* extract-cgs-G.rr: Nabeshima による CGS 計算プログラム kcgs_main の出力から断片とGroebner基底を抽出するスクリプト。
* extract-cgs-g-script.rr: extract-cgs-G.rr を用いて CGS から断片とGroebner基底を抽出するスクリプト。G2-reverse.dat を用いて逆順（
* zero-dimensional-test.rr: CGS の(断片,基底)のペアに収められている基底 G で生成されるイデアル \<G\> がゼロ次元かどうかを判定する。実際には G に属する多項式 g の頭単項式
 (leading monomial) が1変数であり、すべての変数が網羅されているかどうかを調べる。
* generate-G2-new.rr: G2-reverse.dat の CGS に対し、以下の処理を行う。
    1. 断片が実根を含むもののみを抽出。
    1. sqrt(2) の定義多項式　a^2-2 を取り除く。
    1. 第 k 番目の断片のGroebner 基底 G_k が s_1 を頭単項式に持つ多項式を持たない場合は以下の処理を行う。
        1. g = c_1^2 + s_1^2 = 1 を含む場合（c_1 を頭単項式に含む多項式が g の場合）: g を G_k から取り除く。
        1. それ以外の場合（c_1 を頭単項式に含む多項式が上の g でない場合）: g を G_k に加える。
* generate-G2-new-a.rr: generate-G2-new.rr から「sqrt(2) の定義多項式　a^2-2 を取り除く」を除いたもの。
* zero-dimensional-test.rr: CGS（およびその断片）のGroebner基底の各多項式のleading monomialを出力。
* zero-dimensional-test-2.rr: CGS（およびその断片）のGroebner基底の各多項式のleading monomialとleading coefficientを出力。
* zero-dimensional-test-script-G2-new.rr: G2-new.dat のCGSに対してzero-dimensional-test.rrの計算を実行。
* zero-dimensional-test-2-script-G2-new-a.rr: G2-new-a.dat のCGSに対してzero-dimensional-test-2.rrの計算を実行。
* generate-G2-2.rr: G2-new.dat, G2-new-a.dat からそれぞれ G2-2.dat, G2-2-a.dat を生成する。

## データファイル

* G0.dat: 単項式順序 DegRevLex で計算したCGS. (断片,基底)のペアは計算順に収められている。成分のペアは最後に計算されたものが最初に入っているので注意。
* G2.dat: 辞書式順序 Lex で計算したCGS. (断片,基底)のペアは計算順に収められている。成分のペアは最後に計算されたものが最初に入っているので注意。
* G2-reverse.dat: G2.dat の(断片,基底)のペアの並び方を逆にしたもの。成分のペアのは計算した順番に収められている。
* G2-new.dat: G2-reverse.dat から generate-G2-new.rr による処理を行ったCGSのデータ。
* G2-new-a.dat: G2-new.dat と同じ内容だが a (= sqrt(2)) をそのまま残している。
* G2-2.dat: G2-new.dat から generate-G2-2.rr による処理を行ったCGSのデータ。次の Hermite の2次形式の計算に用いる。
* G2-2-a.dat: G2-2.dat と同じ内容だが a (= sqrt(2)) をそのまま残している。