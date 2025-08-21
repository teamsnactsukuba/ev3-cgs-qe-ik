# algn

sqrt(2) を代数的数として扱った計算を記録

## スクリプト

* cgs-F-input-algn.0.rr: 単項式順序 DegRevLex でCGSを計算するスクリプト。出力先はG0.dat。
* cgs-F-input-algn.2.rr: 単項式順序 Lex でCGSを計算するスクリプト。出力先はG2.dat。
* extract-cgs-G.rr: Nabeshima による CGS 計算プログラム kcgs_main の出力から断片とGroebner基底を抽出するスクリプト。
* extract-cgs-g-script.rr: extract-cgs-G.rr を用いて CGS から断片とGroebner基底を抽出するスクリプト。G2-reverse.dat を用いて逆順（
* zero-dimensional-test.rr: CGS の(断片,基底)のペアに収められている基底 G で生成されるイデアル \<G\> がゼロ次元かどうかを判定する。実際には G に属する多項式 g の頭単項式
 (leading monomial) が1変数であり、すべての変数が網羅されているかどうかを調べる。
* arrange-G.rr: G2-reverse.dat の CGS に対し、以下の処理を行う。
    1. 断片が実根を含むもののみを抽出。
    1. sqrt(2) の定義多項式　a^2-2 を取り除く。
    1. 第 k 番目の断片のGroebner 基底 G_k が s_1 を頭単項式に持つ多項式を持たない場合は以下の処理を行う。
        1. g = c_1^2 + s_1^2 = 1 を含む場合（c_1 を頭単項式に含む多項式が g の場合）: g を G_k から取り除く。
        1. それ以外の場合（c_1 を頭単項式に含む多項式が上の g でない場合）: g を G_k に加える。
* arrange-G-a.rr: arrange-G.rr から「sqrt(2) の定義多項式　a^2-2 を取り除く」を除いたもの。

## データファイル

* G0.dat: 単項式順序 DegRevLex で計算したCGS. (断片,基底)のペアは計算順に収められている。成分のペアは最後に計算されたものが最初に入っているので注意。
* G2.dat: 辞書式順序 Lex で計算したCGS. (断片,基底)のペアは計算順に収められている。成分のペアは最後に計算されたものが最初に入っているので注意。
* G2-reverse.dat: G2.dat の(断片,基底)のペアの並び方を逆にしたもの。成分のペアのは計算した順番に収められている。
* G2-new.dat: G2-reverse.dat から arrange-G.rr による処理を行ったCGSのデータ。次の Hermite の2次形式の計算に用いる。
* G2-new-a.dat: G2-new.dat と同じ内容だが a (= sqrt(2)) をそのまま残している。