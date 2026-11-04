# Deterministic basis に対する行列表現

`src.arithmetic.torsion.ChangeOfBasis(B1, B2)` は、SQIsign-v2 と同じ行の規約で
`B2 = C * B1` を満たす行列を返します。点対は列ベクトルと解釈します。

```text
C = [[a, b], [c, d]]
B1 = (P, Q)
B2 = (a*P+b*Q, c*P+d*Q)
```

係数環は `Zmod(2^f)`、現在は `f=240` です。秘密鍵の奇数位数の混合行列
`M` は従来どおり列に像を格納し、法も `N` です。両者を転置なしで混用しません。

## 鍵・証明の形式

```python
pk = (E, R, S, RmS, hint)
sk = (E, I, M, C0, U_can, V_can)
pi = (E3, C3, hint3)
```

`B0` は公開前計算の `(U0_Fp4,V0_Fp4)` を固定した同型写像で twist に移した
決定的な基底です。`standard_basis()` が返します。
`Bcan(E,hint)` と `Bcan(E3,hint3)` は Montgomery の基底生成と符号規約で
再現する完全な `2^f`-torsion 基底です。`A=0` では hint を使わず、決定的な点探索を使います。

旧形式の点対との関係は次です。

```text
旧秘密鍵の (a*U0+b*V0, c*U0+d*V0) = C0 * B0
秘密鍵の   (U_can,V_can)           = Bcan(E,hint) （点対をそのまま保持）
旧証明の   (U3,V3)                 = C3 * Bcan(E3,hint3)
```

Keygen は `C0 = ChangeOfBasis((U,V), Bcan(E,hint))` を計算します。
ここで `(U,V)` は、従来どおり秘密イデアルのノルムの逆数で正規化した同種写像の像です。
canonical 基底 `(U_can,V_can)` は点対として秘密鍵に保存し、Evalで直接使います。
この基底用の変換行列と、Evalでの基底の再生成は不要です。

Eval は `C0` を `IdealToIsogenyFromMatrix` に渡し、標準基底の像へ直接行列を
作用させます。`C0*B0` の点への展開と、それを元に戻す2回の PointCoords を省きます。
基底のペアリングも `e(C0*B0)=e(B0)^det(C0)` から確認します。
証明の生成時にだけ `ChangeOfBasis(Bcan(E3,hint3),(U3,V3))` を計算します。
秘密鍵のpublic側は canonical 基底そのものなので、証明の行列への追加の基底補正は不要です。

Verify は `C3` を点へ展開し、従来と同じ検証を行います。
行列のサイズ・係数環・奇数行列式、hint の範囲、復元基底の位数と独立性を確認し、
不正な証明は `False` にします。行列は完全な `2^240` の点対を保持し、
検証の内部ではこれまでどおり2倍して `2^239`-torsion に縮小します。

```python
from src.encoding.encode import decode_secret_key, decode_proof

E, I, M, U0, V0, U_can, V_can = decode_secret_key(key.pk, key.sk)
E3, U3, V3 = decode_proof(evaluation.pi)
```

通常・償却評価の公開 API の呼び方は変わりません。評価キャッシュは、不要な差分点を除いた
`(K,Fx)` です。旧7要素の秘密鍵、点対を含む旧証明、旧評価キャッシュは入力形式が変わります。
公開鍵も旧4要素から5要素に変わり、hint は `pk[4]` に置かれます。
これは Sage オブジェクトの表現であり、バイト列へのシリアライズは実装していません。
検証APIは `Verify(pk,pi,x,v)` です。公開鍵検査と基底準備を毎回行い、
検証用キャッシュや償却検証のAPIは公開しません。

## Keygen の IdealToIsogeny の返り値

`src.protocol.keygen` の `IdealToIsogeny` は backend の `IdealToIsogenyKeygen` の別名です。
返り値は次の8要素です。

```python
success, curve, P, Q, PmQ, U, V, UmV = IdealToIsogeny(I_sk)
```

`P,Q,PmQ` は非smoothな `N`-torsion の像、`U,V,UmV` は `2^f`-torsion の像です。
6点はいずれも `KummerLine(curve)` 上の射影座標 `(X:Z)` を持つ `KummerPoint` です。
`PmQ` と `UmV` は、対応する点対の差のx座標を表します。
奇数位数の像も `lift=False` の同種写像の出力から、正規化せずに受け取ります。

Keygen は `P,Q,PmQ` と混合行列から、差分演算で `x(R),x(S),x(R-S)` を計算します。
係数は `N` を法として処理します。メモリ上の公開鍵はKummerの差分3点を保持し、
バイト列には `x(R),x(S)` と差分点を選ぶ1バイトのみ保存します。
デコード時に双二次関係から平方根を1回計算して `x(R-S)` を復元します。
Eval 用の `IdealToIsogeny` / `IdealToIsogenyFromMatrix` の返り値は
従来どおり `(success,curve,U,V,UmV)` です。

## Keygen から Verify までの射影演算

公開鍵の `R,S,RmS` はx座標の射影表現です。Eval・Verifyは
`odd_linear_combination(R,S,RmS,r,s,N)` で入力核を計算します。
公開鍵の検査も `(X:Z)` のまま行い、点の所属、非零の `N`-torsion、
差分3点の双二次関係を確認します。`(0:0)` と不正な差分点は拒否します。
この検査自体はWeilペアリングによる独立性検査を行いません。
生成された基底の独立性とペアリング規約は回帰テストで別途照合します。

符号付きの完全な点を使う経路は `src.arithmetic.projective.JacobianPoint` へ移しました。
`y^2=x^3+A*x^2+B*x+C` に対して `x_short=x+A/3` とし、内部では短い
Weierstrassモデルの `(X:Y:Z)`、すなわち `x_short=X/Z^2, y=Y/Z^3` を保持します。
加算・倍算・符号付きスカラー倍の演算中には除算や平方根を使いません。
長いスカラーには幅4の符号付き窓を使います。

このアダプタを通すのは、IdealToIsogenyの自己準同型行列、基底変換行列、
表現の `N` 倍とcofactor倍、積同種写像のkernelの倍算・加算です。
PointCoordsの再構成チェックもJacobian演算です。
smooth基底の位数と独立性は、Kummer上で2-torsionまで倍算し、
2点が異なる非零の2-torsionになることから確認します。

残る正規化・点の復元は次の境界にあります。

- WeilペアリングとPARIのPointCoordsへ渡す完全な点。
- thetaのgluing行列が要求する正規化済みの座標。
- Kummerの入力核・参照点を符号付きtheta APIに渡す際のlift。
- 秘密鍵に保存するcanonical点対、および行列へ符号化するsmooth像のy座標復元。

したがって全経路で有限体の除算や平方根がなくなったわけではありません。
上流のソースは変更せず、ローカルの算術アダプタで射影演算へ移しています。
小さい有限体上の全加算対、長いスカラー、射影再スケーリングを
`tests/test_projective.py` でSageの楕円曲線演算と照合します。

### thetaコードが受け取る座標

`backend.CouplePoint` は、上流の `CouplePoint` に `JacobianPoint` の
オブジェクトを格納します。thetaのkernel処理が呼ぶ `P+Q` や `n*P` は
Jacobian演算になります。一方、上流のgluingコードは `T[0]`、`T[2]` を読み、
非零点では `T[2]==1` を仮定しています。

このため `JacobianPoint.__getitem__` は内部の `(X,Y,Z)` を返さず、

```text
x = X/Z^2 - A/3
y = Y/Z^3
(x, y, 1)
```

という正規化済みの元の曲線上の座標を返します。零点はSageと同じ `(0,1,0)` です。
内部のJacobian座標をそのままgluingの入力座標として解釈させると、
`x=X/Z` という上流の前提と一致しません。このアダプタは座標の読み出し時に変換し、
上流の式と非零点の `T[2]==1` の検査を成立させます。

gluingを通った後の点は上流の `ThetaPoint` で、4個のtheta座標を使います。
したがって、Jacobian座標を使うのは符号付き楕円曲線点の演算部分です。
theta同種写像全体をJacobian座標で計算する実装ではありません。

## bi_dlp から PointCoords への置き換え

以前の `a,b,c,d = bi_dlp(P,Q,R,S,n)` は、次の呼び出しで置き換えられます。

```python
from src.arithmetic.torsion import PointCoords, point_coords_pairing

zeta = point_coords_pairing(R, S, order=n)
a, b = PointCoords(P, R, S, zeta, order=n)
c, d = PointCoords(Q, R, S, zeta, order=n)
# P = a*R+b*S, Q = c*R+d*S  (mod n)
```

`order=n` は、既に位数 `n` の基底 `(R,S)` を使い、追加のスカラー倍を行わない指定です。
既存の `e=e` は、位数 `2^f` の基底から `2^(f-e)` 倍して座標を求める指定です。
両者は同時には指定できません。省略時の `e=f` の動作は同じです。
`ChangeOfBasis(B1,B2,order=n)` は、この規約の行列を `Zmod(n)` 上で返します。
`order` を使う場合は `B1` を明示します。

`precompute.endomorphism.action_matrix` も上記のPointCoordsを使います。
前計算の行列は従来どおり像を列に格納し、上流の行の規約へは `params.py` で転置します。
旧 `src.arithmetic.torsion.bi_dlp` と `precompute.endomorphism.bi_dlp` は、
同じPointCoordsを呼ぶ互換入口として保持します。Sageの汎用 `discrete_log` による
重複実装は削除しました。

入力はy座標を持つSage点または `JacobianPoint` です。PARIへ渡す境界で正規化し、
座標による点の再構成チェックはJacobian演算で行います。
Kummer点を使う場合は、差分点とともに点対を復元してから渡します。

## 上流 ec.PointCoords をそのまま呼べない理由

「PointCoords 自体が使えない」わけではありません。
ローカルの `src.arithmetic.torsion.PointCoords` は、上流と同じ座標・スケーリング規約で使えます。
相違はペアリングの計算方法です。

上流 `ec.PointCoords` は `pari.elltatepairing` の値を
`(p^2-1)/2^e` 乗し、2本の離散対数を座標として返します。
PARI が返す値は **非簡約 Tate ペアリング**で、体 `k` の
`k*/(k*)^n` の類です。一方、PARI の Weil ペアリングは交代双線形写像です。
([PARI公式仕様](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html#elltatepairing),
[Weilの仕様](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html#ellweilpairing))

上流の標準モデルは二次体 `Fp2` 上ですが、本実装は奇数位数 `N` の点も扱うため
twist の `Fp4` モデルを使います。`Fp2` 用の最終指数をそのまま適用しても、
`Fp4` の Tate の類から必要な双線形写像は得られません。
今回の再現では、これが単なる想定の違いではなく実際の係数の誤りになることを確認しました。

`n=2^240`、`h=2^239`、固定標準基底を `(P,Q)` とし、**上流コードの指数で計算した値**を
`T` と書くと、実測は次です。

```text
zeta = T(P,Q) は位数 n
T(P,P) = T(Q,Q) = -1
T(Q,P) = zeta^(-1)
T(P+Q,Q) != T(P,Q)*T(Q,Q)
```

最後の式のとおり、この `T` は双線形性を満たしません。
自己ペアリングの項だけを2×2行列で補正すれば一般に直る、という状況でもありません。
最小の失敗例は `X=P` です。上流版の2本目の離散対数が `h` になり、
`(1,0)` の代わりに `(1,h)` を返そうとして最後の再構成 `assert` に失敗します。
一方、`X=3P+5Q` では正しく `(3,5)` が得られ、少数の成功例だけでは問題を検出できません。

旧秘密鍵の点でも最上位ビットのずれを観測しました。
旧診断テストでは乱数状態による入力の変化を避け、固定入力 `X=P` で
Weil版の `(1,0)` と上流Tate版の `(1,h)` を比較しました。
本体が使わない上流Tate版の失敗再現テストは削除し、この診断結果のみ記録しています。

### 指数の変更だけで解決するか

四次体用の `(p^4-1)/n` に変えると、この標準基底のペアリングの位数は
`n` から **`n/2`** に落ちます。このことも旧診断テストで確認しました。
元の `zeta` を使えば、新しい値は `zeta^(p^2+1)` です。
`p` は奇数なので `v2(p^2+1)=1` となり、最上位の1ビットを区別できなくなります。
したがって、最終指数だけ変えて元の `fflog(...,n)` に渡す方法では完全な座標を復元できません。

一方、点を固定同型写像で `E0` に移し、二次部分体上の点・曲線として作り直すと、
同じ失敗入力でも上流 PointCoords が正しい係数を返します。
つまり Tate 版を利用する余地はありますが、曲線と点のモデル・体の輸送まで含む
アダプタが必要です。現状の各曲線・基底へそのまま置換することはできません。

### 現在の実装

`src.arithmetic.torsion.PointCoords` は同梱の
`theta.utilities.discrete_log.BiDLP` を再利用し、PARI の Weil ペアリングと
`fflog` で座標を求めます。完全な基底なら

```text
zeta = e_n(P,Q)
e_n(aP+bQ,Q)  = zeta^a
e_n(aP+bQ,-P) = zeta^b
```

が成り立ち、`Fp4` 上でも法 `n` の座標を240ビットすべて復元できます。
`zeta` 引数を共有する場合は `point_coords_pairing` の Weil 値を使います。
上流の Tate 値とは互換ではありません。
`e<f` の API は `X=2^(f-e)*(aP+bQ)` を返す規約です。
ローカル ChangeOfBasis もこの規約に従います。

本プロトコルの座標計算の対象は滑らかな **2べき位数**です。
`PointCoords` の引数 `e` は `2^e` の指数です。`order=n` は小さい奇数・合成数位数にも使えますが、
239ビット素数位数 `N` の一般の離散対数を高速に解く機能ではありません。
奇数位数の自己準同型作用は、従来どおり Frobenius 固有基底の既知の関係から計算します。

回帰テストは次で実行します。

```sh
sage -python -m tests.test_protocol
sage -python -m tests.test_pointcoords
```
