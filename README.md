# VarQEC

This is the code repository for paper "Quantum variational learning for quantum error-correcting codes" [arxiv-link](https://arxiv.org/abs/2204.03560)

quickstart

```bash
python draft00.py
```

You could chnnge `qecc_str` in `draft00.py` to find various QECC.

set `qecc_str=((5,2,de(2)=3))[5]` to find asymmetrical QECC.

The code we found and the corresponding random seed

| code | status | seed | note |
| :-: | :-: | :-: | :-: |
| `((5,2,3))[5]` | pass | 2585155996 | non-degenerate, about `50` step |
| `((5,2,de(2)=3))[5]` | pass | 391718720280628643 | about `80` step |
| `((5,6,2))[5]` | pass | - | |
| `((5,6,2))[6]` | pass | 6454016727 | |
| `((6,2,3))[5]` | pass | 9989362660 | |
| `((6,2,de(2)=4))[5]` | - | - | |
| `((6,4,de(2)=3))[3]` | pass | 5498500481 | |
| `((6,2,de(0.5)=2))[5]` | pass | 263168234357481788 | |
| `((7,2,3))[5]` | pass | 421928941517207614 | |
| `((7,16,2))[5]` | pass | 305656324415284766 | |
| `((7,3,de(0.5)=2))[5]` | pass | 861426403674806629 | |
| `((7,8,de(2)=3))[3]` | pass | 496387136 | |
| `((8,8,3))[5]` | pass | 63624525 | |
| `((8,3,de(2)=4))[3]` | - | - | |
| `((9,8,3))[5]` | pass | 4834842764 | |
| `((11,32,3))[5]` | pass | 8230370079 | |
| `((12,16,3))[5]` | pass | 6102247691 | |
| `((12,64,3))[5]` | pass | 626895936804222064 | |
| `((13,128,3))[5]` | pass | 435676291134251631 | |
| `((14,256,3))[5]` | pass | 395281354919490522 | |

above `((n,K,d))[L]`

Some of `seed` might be not correct but that's time-consuming to verify, please report if you find some seed is wrong.
