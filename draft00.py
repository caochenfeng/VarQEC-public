import numpy as np
import scipy.optimize
import qiskit

from utils import make_error_list, make_asymmetric_error_set, hf_callback_wrapper
from utils import parse_qecc_str, make_bipartite_connection

Aer_simulator = qiskit.Aer.get_backend('statevector_simulator')


def run_circuit(num_qubit, num_layer, params, code_ind, connection_list):
    qr = qiskit.QuantumRegister(num_qubit)
    circ = qiskit.QuantumCircuit(qr)

    tmp0 = [x for x,y in enumerate(bin(code_ind)[-1:1:-1]) if y=='1']
    for x in tmp0:
        circ.x(x)

    params_s = params[:2 * num_qubit * (num_layer + 1)].reshape([num_layer + 1, num_qubit, 2])
    params_d = params[2 * num_qubit * (num_layer + 1):].reshape([num_layer, len(connection_list)])
    for ind0 in range(num_layer):
        for ind1 in range(num_qubit):
            circ.rx(params_s[ind0,ind1,0], qr[ind1])
            circ.rz(params_s[ind0,ind1,1], qr[ind1])
        for ind1 in range(len(connection_list)):
            circ.rzz(params_d[ind0,ind1], qr[connection_list[ind1,0]], qr[connection_list[ind1,1]])
    for ind1 in range(num_qubit):
        circ.rx(params_s[-1,ind1,0], qr[ind1])
        circ.rz(params_s[-1,ind1,1], qr[ind1])
    circ = qiskit.transpile(circ, Aer_simulator)
    ret = Aer_simulator.run(circ).result().data(0)["statevector"]
    return ret


def loss_function_wrapper(num_qubit, num_logical_dim, error_list, num_layer, connection_list, loss_type='L2', sample_ratio=None, seed=None):
    assert loss_type in {'L1', 'L2'}
    K = num_logical_dim
    np_rng = np.random.default_rng(seed)
    def hf_loss(params):
        code = [run_circuit(num_qubit, num_layer, params, x, connection_list) for x in range(K)]
        if sample_ratio is None:
            sub_error_list = error_list
        else:
            N0 = len(error_list)
            tmp0 = np_rng.choice(N0, min(N0, int(N0*sample_ratio)), replace=False, shuffle=False)
            sub_error_list = [error_list[x] for x in tmp0]
        loss = 0
        if loss_type=='L2':
            for e in sub_error_list:
                loss += sum([np.abs(np.vdot(code[i], e.dot(code[j])))**2 for i in range(K-1) for j in range(i+1, K)])
                tmp0 = np.array([np.vdot(x, e.dot(x)) for x in code])
                loss += np.sum(np.abs(tmp0 - tmp0.mean())**2)/2
        else:
            for e in sub_error_list:
                loss += sum([np.abs(np.vdot(code[i], e.dot(code[j]))) for i in range(K-1) for j in range(i+1, K)])
                tmp0 = np.array([np.vdot(x, e.dot(x)) for x in code])
                loss += np.sum(np.abs(tmp0 - np.mean(tmp0)))/2
        return loss
    num_parameter = 2 * num_qubit * (num_layer + 1) + num_layer * len(connection_list)
    return hf_loss, num_parameter


qecc_str = '((5,2,3))[5]' #((n,K,d))[L]
# qecc_str = '((5,2,de(2)=3))[5]'

tmp0 = parse_qecc_str(qecc_str)
num_qubit = tmp0['num_qubit']
num_logical_dim = tmp0['num_logical_dim']
distance = tmp0['distance']
num_layer = tmp0['num_layer']
weight_z = tmp0['weight_z']

connection_list = make_bipartite_connection(5, 2)
if weight_z is None:
    error_list = make_error_list(num_qubit, distance)
else:
    error_list = make_asymmetric_error_set(num_qubit, distance, weight_z)

hf_loss_L2, num_parameter = loss_function_wrapper(num_qubit,
        num_logical_dim, error_list, num_layer, connection_list, loss_type='L2')
hf_loss_L1, num_parameter = loss_function_wrapper(num_qubit,
        num_logical_dim, error_list, num_layer, connection_list, loss_type='L1')

seed = 2585155996
np_rng = np.random.default_rng(seed)

theta0 = np_rng.uniform(0, 2*np.pi, num_parameter)
hf_callback = hf_callback_wrapper(hf_loss_L2, print_freq=10)
theta_optim = scipy.optimize.minimize(hf_loss_L2, theta0, method="BFGS", tol=1e-6, callback=hf_callback)

# optimize L1-norm loss function if needed
# theta_optim = scipy.optimize.minimize(hf_loss_L1, theta_optim.x, tol=1e-8)
print(theta_optim.fun, theta_optim.x)
