import time
import itertools
import numpy as np

def hf_split_element(np0, num_of_each):
    # split np0 (len(np0)=N) into len(num_of_each) sets, elements in each set is unordered
    assert len(num_of_each) > 0
    assert sum(num_of_each) <= len(np0)
    tmp0 = list(range(np0.shape[0]))
    tmp1 = np.ones(np0.shape[0], dtype=np.bool_)
    if num_of_each[0]==0:
        if len(num_of_each)==1:
            yield (),
        else:
            tmp0 = (),
            for x in hf_split_element(np0, num_of_each[1:]):
                yield tmp0+x
    else:
        for ind0 in itertools.combinations(tmp0, num_of_each[0]):
            ind0 = list(ind0)
            if len(num_of_each)==1:
                yield tuple(np0[ind0].tolist()),
            else:
                tmp1[ind0] = False
                tmp2 = np0[tmp1]
                tmp1[ind0] = True
                tmp3 = tuple(np0[ind0].tolist()),
                for x in hf_split_element(tmp2, num_of_each[1:]):
                    yield tmp3 + x
# z0 = hf_split_element(np.arange(num_qubit), [0,1,0])

def make_asymmetric_error_set(num_qubit, distance, weight_z=1):
    # 0<=nx+ny+nz<=num_qubit
    # 0<=nx,ny,nz
    # nx+ny+cz*nz<distance
    assert weight_z>0
    I = np.eye(2)
    sx = np.array([[0.0, 1.0], [1.0, 0.0]])
    sy = np.array([[0.0, -1j], [1j, 0.0]])
    sz = np.array([[1.0, 0.0], [0.0, -1.0]])

    ret = []
    for nxy in range(min(num_qubit,distance)):
        tmp0 = int(np.ceil((distance-nxy)/weight_z))
        for nz in range(min(num_qubit-nxy+1, tmp0)):
            if (nxy==0) and (nz==0):
                continue
            for nx in range(nxy+1):
                ny = nxy - nx
                for indx,indy,indz in hf_split_element(np.arange(num_qubit), [nx,ny,nz]):
                    tmp = [I] * num_qubit
                    for x in indx:
                        tmp[x] = sx
                    for x in indy:
                        tmp[x] = sy
                    for x in indz:
                        tmp[x] = sz
                    ret.append(kron_list(tmp))
    return ret


def kron_list(op_list):
    ret = op_list[0]
    for x in op_list[1:]:
        ret = np.kron(ret, x)
    return ret


def make_error_list(n, d):
    I = np.eye(2)
    X = np.array([[0,1], [1,0]])
    Y = np.array([[0,-1j], [1j,0]])
    Z = np.array([[1,0], [0,-1]])
    ret = []
    for wt in range(1, d):
        for c in itertools.combinations(np.arange(n), wt):
            for p in itertools.product([X, Y, Z], repeat=wt):
                tmp = [I] * n
                for x in range(wt):
                    tmp[c[x]] = p[x]
                ret.append(kron_list(tmp))
    return ret


def hf_callback_wrapper(hf_fval, state:dict=None, print_freq:int=1):
    if state is None:
        state = dict()
    state['step'] = 0
    state['time'] = time.time()
    state['fval'] = []
    state['time_history'] = []
    def hf0(theta):
        step = state['step']
        if (print_freq>0) and (step%print_freq==0):
            t0 = state['time']
            t1 = time.time()
            fval = hf_fval(theta)
            print(f'[step={step}][time={t1-t0:.3f} seconds] loss={fval}')
            state['fval'].append(fval)
            state['time'] = t1
            state['time_history'].append(t1-t0)
        state['step'] += 1
    return hf0


def parse_qecc_str(qecc_str):
    assert (qecc_str[:2]=='((') and ('))[' in qecc_str) and (qecc_str[-1]==']')
    tmp0 = qecc_str[2:-1]
    tmp1,tmp2 = tmp0.split('))[',1)
    tmp3 = tmp1.split(',',2)
    num_qubit = int(tmp3[0])
    num_logical_dim = int(tmp3[1])
    if '=' in tmp3[2]:
        # ((6,2,de(2)=4))[8]
        weight_z = float(tmp3[2].split('(',1)[1].split(')',1)[0])
        distance = int(tmp3[2].split('=',1)[1])
    else:
        weight_z = None
        distance = int(tmp3[2])
    num_layer = int(tmp2)
    ret = dict(num_qubit=num_qubit, num_logical_dim=num_logical_dim,
            weight_z=weight_z, distance=distance, num_layer=num_layer)
    return ret


def make_bipartite_connection(num_qubit, num_logical_dim):
    num_logical_qubit = int(np.ceil(np.log2(num_logical_dim)))
    ret = np.array([(x,y) for x in range(num_logical_qubit) for y in range(num_logical_qubit, num_qubit)])
    return ret
