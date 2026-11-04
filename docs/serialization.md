# バイト列形式

`Keygen.pk`, `Keygen.sk`, `Eval.pi`, `Eval.v`（Amortized_Eval も同様）は
`bytes` です。Eval の入口で公開鍵と秘密鍵を、Verify の入口で公開鍵・証明・
出力をデコードします。Amortized_Eval は返した `state` にデコード済み鍵を
保持し、同じ鍵で `state` を再利用する2回目以降は鍵のデコードを省略します。
別の鍵での再利用はエラーです。従来の2要素タプルを手動で作った場合や、
`amortized=(K,Fx)` のみ渡す場合は鍵をデコードします。
Verify は不正なバイト長・範囲・行列などを拒否して
`verify=False` を返します。デコードと基底復元も計測時間に含まれます。
`amortized` / `state` は実行中の算術オブジェクトのままです。
`state` は従来どおり `(I_aux,E_aux)` として添字参照・展開できる
`AmortizedState` です。鍵キャッシュを保持するため、そのオブジェクトを
`state=previous.state` の形で渡してください。

参照仕様は指定された **SQIsign, 2025-07-07, §4.6 Binary format**、
イデアル生成元は Algorithm 3.8 です。これは DeuringVUF の形式であり、
SQIsign の鍵や署名とのバイト単位の相互運用形式ではありません。

## 基本型

- パラメータセットは実行前に共有します。ヘッダやパラメータ ID は含めません。
- `L = ceil(bit_length(p)/8)`、`T = ceil(f/8)`、`D = ceil(bit_length(N)/8)`。
  Fp の幅は §4.6 の最小幅に従います。付録 B の64ビット limb 丸めは使いません。
- Fp は `[0,p)` の固定幅 L バイト unsigned little-endian。Fp² は
  `a + b*i`（`i²=-1`）を `LE(a) || LE(b)` にします。現在の Fp4 内の
  Fp² 埋め込みを `context.zeta4` で変換し、Fp² の外の値は拒否します。
- 曲線は `y²=x³+A*x²+x` の A のみで表します。
- non-smooth basis は `x(R), x(S)` の Fp² 表現と、差分点を選ぶ1バイトの
  `branch`（0 または1）で保存します。無限遠点は許可しません。
  射影座標のスケーリングに依存しません。
- hint は `hA | (h << 1)` の1バイト。A=0 は独自拡張としてゼロを置き、
  デコード後は `None` として決定的基底を生成します。
- 行列は行優先の `a || b || c || d`。剰余環の係数は最小非負剰余の
  unsigned little-endian です。各行が点の係数を表し、可逆性を検査します。
- イデアルの norm と quaternion の整数係数は仕様の Z に従う固定幅 L の
  signed little-endian（2の補数）です。範囲外はエラーにします。

## 連結順

| オブジェクト | 連結する値 | バイト長 |
| --- | --- | --- |
| 公開鍵 | `A || hint || x(R) || x(S) || branch` | `6L+2` |
| 秘密鍵 | `pk || norm(I) || gamma[0..3] || M_N || C0` | `len(pk)+5L+4D+4T` |
| 証明 | `A3 || C3 || hint3` | `2L+4T+1` |
| 出力 v | `j(Ex)` | `2L` |

公開鍵は SQIsign の `A || hint` に圧縮した odd basis を追加しています。
現在のパラメータでは公開鍵194バイト、秘密鍵594バイト、証明185バイト、
出力64バイトです。旧形式の公開鍵・秘密鍵とは互換性がありません。

差分点は y 座標を復元せず、次の二次方程式から復元します。
`r=x(R), s=x(S), d=(r-s)^2, m=(rs+1)(r+s)+2Ars` とすると、
候補は `(m+sqrt(Δ))/d` と `(m-sqrt(Δ))/d`、`Δ=m²-(rs-1)²d` です。
平方根は Fp² 内にあることを検査し、`encode_fp2(root)` と
`encode_fp2(-root)` をバイト列として辞書順比較して小さい側に符号を固定します。
`branch=0` は加算側、`branch=1` は減算側を選びます。
`r=s`、重根、Fp² に平方根がない入力、0・1以外の branch は拒否します。
復元後の API は従来どおり `(E,R,S,RmS,hint)` を返します。

秘密鍵は SQIsign の公開鍵・イデアル・行列の構成に、DeuringVUF の `M_N`
（係数 D バイト）を加えています。C0/C3 の係数は T バイトです。
証明は SQIsign の署名ではなく、現在の `(E3,C3,hint3)` を直接表現します。

イデアルは `I = O0*gamma + O0*norm(I)` として復元します。Algorithm 3.8 で
生成元を選び、norm と互いに素な分母を払います（§4.6 脚注1）。整数係数を
norm の法で中央剰余にして幅を抑え、エンコード時に元のイデアルと一致する
ことを確認します。デコード時にも norm を検査します。
秘密鍵中の `U_can,V_can` は格納せず、埋め込んだ公開鍵の hint から復元します。
公開鍵が別途渡された場合は秘密鍵内の公開鍵と一致することも確認します。

## API

```python
from src.encoding import serialization as wire

pk = wire.decode_public_key(key.pk)          # (E,R,S,RmS,hint)
sk = wire.decode_secret_key(key.sk, key.pk)  # (E,I,M_N,C0,U_can,V_can)
pi = wire.decode_proof(evaluation.pi)        # (E3,C3,hint3)
j = wire.decode_output(evaluation.v)        # Fp4 内の Fp² の元

assert wire.encode_public_key(pk) == key.pk
assert wire.encode_proof(pi) == evaluation.pi
```

各型に対応する `encode_*` / `decode_*` は `src/encoding/serialization.py` にあります。
`encode_secret_key(pk,sk)` はデコード済みの2つのタプルを受け取ります。
既存の `src.encoding.encode.encode_proof/decode_proof` は引き続き点対と行列の変換用で、
バイト列用関数とは別です。互換性のため Eval/Verify は従来の算術オブジェクトも
受け付けます。バイト列の decode は構文と局所的な条件を検査する処理であり、
公開鍵の torsion・差分関係や証明の有効性の検査は Verify が行います。

テスト:

```sh
sage -python -m tests.test_serialization
```
