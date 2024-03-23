# Import necessary libraries and modules
from qiskit import QuantumCircuit, transpile, assemble, IBMQ, QuantumRegister, ClassicalRegister, execute, Aer
from qiskit.circuit.library import QFT
from qiskit_ibm_runtime import QiskitRuntimeService

from Crypto.PublicKey import RSA
import binascii
from typing import Optional, Dict, Any, List, Counter
from math import gcd, pi
from numpy.random import randint
import numpy as np
import math
import random
import time

# Get the least busy quantum device
# This is done by filtering the available backends based on the number of qubits, whether it's a simulator, and its operational status
provider = QiskitRuntimeService()
backends = provider.backends(filters=lambda x: x.configuration().n_qubits >= 5 and not x.configuration().simulator and x.status().operational==True)
for device in sorted(backends, key=lambda x: x.status().pending_jobs):
    if device.status().operational:
        # Use the device
        break
    else:
        print(f"The device {device.name()} is currently not operational. Trying another device.")

print("Coupling Map Test:", device.configuration().coupling_map)

# Initialize a counter
# This will be used to keep track of incorrect factorizations
counter = 0

# This is the number of times the quantum circuit will be run. More shots will give more accurate results, but will also take longer to run.
shots = 100

'''
NOTE: below is what the encryption/decryption functions *should ideally* look like. When using Shor’s algorithm to factorize the RSA modulus, we need a quantum computer with at least twice as many qubits as the size of the RSA key. In this case, we're using IBM's Eagle R3 quantum processor which is only 127 qubits. So, our RSA key can only realistically be >63 bits for a chance at success. Of course, a key size this low opens us up to chosen plaintext attacks. However, please note this program is only an educational test of Shor's algorith with a QPU. 

from Crypto.Cipher import PKCS1_OAEP
def rsa_encryption(message):
    keyPair = RSA.generate(bits=2048)
    pubKey = keyPair.publickey()
    pubKeyPEM = pubKey.exportKey()
    print(f"Public key:  (n={hex(pubKey.n)}, e={hex(pubKey.e)})")
    print(f"Public key (PEM format):\n{pubKeyPEM}")
    encryptor = PKCS1_OAEP.new(pubKey)
    encrypted = encryptor.encrypt(message)
    print("Encrypted:", binascii.hexlify(encrypted))
    return encrypted, keyPair

def rsa_decryption_with_d(encrypted, d, integer_N):
    # Construct the private key
    private_key = RSA.construct((integer_N, d))
    decryptor = PKCS1_OAEP.new(private_key)
    decrypted = decryptor.decrypt(encrypted)
    return decrypted
'''

def rsa_encryption(message):
    # Generate a 24-bit RSA key pair
    keyPair = RSA.generate(bits=24)
    # Extract the public key from the key pair
    pubKey = keyPair.publickey()
    # Export the public key in PEM format
    pubKeyPEM = pubKey.exportKey()
    print(f"Public key:  (n={hex(pubKey.n)}, e={hex(pubKey.e)})")
    print(f"Public key (PEM format):\n{pubKeyPEM}")
    # Convert the message to an integer
    m = int.from_bytes(message, byteorder='big')
    # Encrypt the message using the public key
    encrypted = pow(m, pubKey.e, pubKey.n)
    print("Encrypted:", binascii.hexlify(encrypted.to_bytes((encrypted.bit_length() + 7) // 8, byteorder='big')))
    return encrypted, keyPair
    
def rsa_decryption_with_d(encrypted, d, integer_N):
    # Construct the private key using the decryption exponent 'd' and the modulus 'integer_N'
    private_key = RSA.construct((integer_N, d))
    # Decrypt the message using the private key
    decrypted = private_key.decrypt(encrypted)
    return decrypted
	
# Encrypt the message
message = b'499'
encrypted, keyPair = rsa_encryption(message)

# Get integer_N and e from the key pair
integer_N = keyPair.n
e = keyPair.e

# Choose a random integer_a that is coprime with integer_N and is in the set {2, 7, 8, 11, 13}
while True:
    integer_a = random.choice([2, 7, 8, 11, 13])
    if math.gcd(integer_a, integer_N) == 1:
        break

def print_circuit_info(circuit):
    print(f"The circuit uses {circuit.num_qubits} qubits.")
    print(f"The circuit has {circuit.size()} operations.")

def aqft(qubits):
    """Performs the Approximate Quantum Fourier Transform on the provided qubits."""
    circuit = QuantumCircuit(len(qubits))
    for j in range(len(qubits)):
        circuit.h(qubits[j])
        for k in range(j+1, min(j+3, len(qubits))):  # Apply controlled rotations only for the next two qubits
            circuit.cp(pi/float(2**(k-j)), qubits[k], qubits[j])
        circuit.barrier()
    print("aqft executed successfully.")
    return circuit

def aiqft(qubits):
    """Performs the Approximate Inverse Quantum Fourier Transform on the provided qubits."""
    # Start with a blank circuit
    circuit = QuantumCircuit(len(qubits))
    # First add swap gates to reverse the order of the qubits
    for i in range(len(qubits)//2):
        circuit.swap(qubits[i], qubits[len(qubits)-i-1])
    # Then apply the AQFT
    for j in range(len(qubits)):
        for m in range(j+1, min(j+3, len(qubits))):
            circuit.cp(-pi/float(2**(m-j)), qubits[m], qubits[j])
        circuit.h(qubits[j])
    print("aiqft executed successfully.")
    return circuit
		
def shors_algorithm(integer_N: int, integer_a: int) -> QuantumCircuit:
    # validate the inputs
    if integer_N < 1 or integer_N % 2 == 0:
        raise ValueError("The input N needs to be an odd integer greater than 1.")
    if integer_a >= integer_N or math.gcd(integer_a, integer_N) != 1:
        raise ValueError('The integer "a" needs to satisfy 1 < a < N and gcd(a, N) = 1.')

    # calculate number of qubits needed
    n = int(np.ceil(np.log2(float(integer_N))))
    m = n

    counting_qubits = QuantumRegister(n, 'counting')
    aux_qubits = QuantumRegister(m, 'auxiliary')
    classical_qubits = ClassicalRegister(n, 'classical')

    shors_circuit = QuantumCircuit(counting_qubits, aux_qubits, classical_qubits)

    # Initialize counting and aux qubits
    shors_circuit.h(counting_qubits)

    # Apply modular exponentiation
    mod_exp_circuit = modular_exponentiation_amod15(counting_qubits, aux_qubits, integer_a)
    shors_circuit.append(mod_exp_circuit, counting_qubits[:] + aux_qubits[:])

    # Print circuit info
    print_circuit_info(shors_circuit)
    
    # Apply AQFT
    shors_circuit.append(aqft(counting_qubits), counting_qubits[:])

    # Apply inverse AQFT
    shors_circuit.append(aiqft(counting_qubits), counting_qubits[:])

    print("shors_algorithm executed successfully.")
    return shors_circuit
	
def run_shors_algorithm(circuit: QuantumCircuit, device, shots: Optional[int] = 100,) -> Dict[str, Any]:
    # Execute the quantum circuit on the specified device
    job = execute(circuit, device, shots=shots)

    # Get the results from the job execution
    result = job.result()

    # Prepare the output dictionary
    out = {
        "measurements": result.get_counts(),
    }

    print("run_shors_algorithm executed successfully.")
    return out
	
def inverse_qft_noswaps(qubits) -> QuantumCircuit:
    """Performs the Inverse Quantum Fourier Transform without swaps on the provided qubits."""
    # Instantiate circuit object
    qft_circuit = QuantumCircuit(qubits)

    # Get number of qubits
    num_qubits = len(qubits)

    # Start on the last qubit and work to the first.
    for k in reversed(range(num_qubits)):
        # Apply the controlled rotations, with weights (angles) defined by the distance to the
        # control qubit. These angles are the negative of the angle used in the QFT.
        # Start on the last qubit and iterate until the qubit after k.
        # When num_qubits==1, this loop does not run.
        for j in reversed(range(1, num_qubits - k)):
            angle = -2 * math.pi / (2 ** (j + 1))
            qft_circuit.cp(angle, qubits[k + j], qubits[k])

        # Then add a Hadamard gate
        qft_circuit.h(qubits[k])

    print("inverse_qft_noswaps executed successfully.")
    return qft_circuit
	
def c_amod15(a: int, power: int) -> QuantumCircuit:
    """Controlled multiplication by a mod 15"""
    if a not in [2,7,8,11,13]:
        raise ValueError("'a' must be 2,7,8,11 or 13")
    U = QuantumCircuit(4)        
    for iteration in range(power):
        if a in [2,13]:
            U.swap(0,1)
            U.swap(1,2)
            U.swap(2,3)
        if a in [7,8]:
            U.swap(2,3)
            U.swap(1,2)
            U.swap(0,1)
        if a == 11:
            U.swap(1,3)
            U.swap(0,2)
        if a in [7,11,13]:
            for q in range(4):
                U.x(q)
    U = U.to_gate()
    U.name = "%i^%i mod 15" % (a, power)
    c_U = U.control()
    end_time = time.time()  # End the timer
    print("c_amod15 executed successfully. ")
    
    return c_U
	
def modular_exponentiation_amod15(counting_qubits, aux_qubits, integer_a) -> QuantumCircuit:
    # Instantiate circuit object
    shors_circuit = QuantumCircuit(counting_qubits, aux_qubits)

    for x in range(len(counting_qubits)):
        r = 2**x
        mod_exp_circuit = c_amod15(integer_a, r)
        shors_circuit.append(mod_exp_circuit, [counting_qubits[x]] + list(aux_qubits[:4]))
        
        # Print circuit info after each call to c_amod15
        print_circuit_info(shors_circuit)
    print("modular_exponentiation_amod15 executed successfully.")
	
def get_factors_from_results(results: Dict[str, Any], integer_N: int, integer_a: int, verbose: bool = True,) -> Dict[str, Any]:
    # Unpack the results dictionary to get the measurement counts
    measurement_counts = results["measurements"]

    # Get the phases from the measurement counts
    phases_decimal = _get_phases(measurement_counts)

    r_guesses = []
    factors = []
    if verbose:
        print(f"Number of Measured phases (s/r) : {len(phases_decimal)}")
    for phase in phases_decimal:
        if verbose:
            print(f"\nFor phase {phase} :")
        # Get the denominator of the phase as an estimate for r
        r = (Fraction(phase).limit_denominator(integer_N)).denominator
        r_guesses.append(r)
        if verbose:
            print(f"Estimate for r is : {r}")
        # Calculate the factors of N using the estimate for r
        factor = [
            math.gcd(integer_a ** (r // 2) - 1, integer_N),
            math.gcd(integer_a ** (r // 2) + 1, integer_N),
        ]
        factors.append(factor[0])
        factors.append(factor[1])
        if verbose:
            print(f"Factors are : {factor[0]} and {factor[1]}")
    # Remove trivial factors (1 and N)
    factors_set = set(factors)
    factors_set.discard(1)
    factors_set.discard(integer_N)
    if verbose:
        print(f"\n\nNon-trivial factors found are : {factors_set}")

    aggregate_results = {"guessed_factors": factors_set}
    
    print("get_factors_from_results executed successfully.")
    return aggregate_results
	
def _get_phases(measurement_counts) -> List[float]:
    # Aggregate the results (i.e., ignore/trace out the query register qubits):
    if not measurement_counts:
        return None

    # First get bitstrings with corresponding counts for counting qubits only (top half)
    num_counting_qubits = int(len(list(measurement_counts.keys())[0]) / 2)

    bitstrings_precision_register = [key[:num_counting_qubits] for key in measurement_counts.keys()]

    # Then keep only the unique strings
    bitstrings_precision_register_set = set(bitstrings_precision_register)
    # Cast as a list for later use
    bitstrings_precision_register_list = list(bitstrings_precision_register_set)

    # Now create a new dict to collect measurement results on the precision_qubits. Keys are given
    # by the measurement count substrings on the register qubits. Initialize the counts to zero.
    precision_results_dict = {key: 0 for key in bitstrings_precision_register_list}

    # Loop over all measurement outcomes
    for key in measurement_counts.keys():
        # Save the measurement count for this outcome
        counts = measurement_counts[key]
        # Generate the corresponding shortened key (supported only on the precision_qubits register)
        count_key = key[:num_counting_qubits]
        # Add these measurement counts to the corresponding key in our new dict
        precision_results_dict[count_key] += counts

    # Convert the bitstrings to decimal to get the phases
    phases_decimal = [_binary_to_decimal(item) for item in precision_results_dict.keys()]

    print("_get_phases executed successfully.")
    return phases_decimal

def _binary_to_decimal(binary: str) -> float:
    fracDecimal = 0
    # Convert fractional part of binary to decimal equivalent
    twos = 2

    for ii in range(len(binary)):
        fracDecimal += (ord(binary[ii]) - ord("0")) / twos
        twos *= 2.0

    # return fractional part
    print("_binary_to_decimal executed successfully.")
    return fracDecimal

def extended_gcd(a, b):
    # This function implements the extended Euclidean algorithm
    # It returns the greatest common divisor of a and b, along with the coefficients of Bézout's identity
    if a == 0:
        return b, 0, 1
    else:
        gcd, x, y = extended_gcd(b % a, a)
        return gcd, y - (b // a) * x, x
    print("extended_gcd executed successfully.")

def get_d(e, phi):
    # This function calculates the multiplicative inverse of e modulo phi
    # This is used in RSA decryption
    gcd, x, _ = extended_gcd(e, phi)
    if gcd != 1:
        raise ValueError("Multiplicative inverse does not exist")
    else:
        return x % phi
    print("get_d executed successfully.")

def run_shors_algorithm_until_success(integer_N: int, integer_a: int, device, shots: Optional[int] = 100, max_attempts: Optional[int] = 20) -> Dict[str, Any]:
    start_time = time.time()  # Start the timer
    for _ in range(max_attempts):
        # Run Shor's algorithm
        shors_circuit = shors_algorithm(integer_N, integer_a)
        print("Circuit:", shors_circuit)
        results = run_shors_algorithm(shors_circuit, device, shots)
        # Get the factors from the results
        aggregate_results = get_factors_from_results(results, integer_N, integer_a, verbose=True)
        factors = aggregate_results["guessed_factors"]

        # If we found the correct factors, return them
        if len(factors) == 2:
            end_time = time.time()  # End the timer
            print(f"Time taken: {(end_time - start_time) * 1000} milliseconds")
            return factors
            print("run_shors_algorithm_until_success executed successfully.")

    # If we didn't find the correct factors after the maximum number of attempts, raise an error
    raise ValueError("Failed to find factors after maximum number of attempts")

# Run Shor's algorithm and get the factors
try:
    factors = run_shors_algorithm_until_success(integer_N, integer_a, device, shots)

    # Check if we found the correct factors
    if len(factors) == 2:
        p, q = factors

        # Compute phi(integer_N)
        phi = (p - 1) * (q - 1)

        # Compute d
        d = get_d(e, phi)

        # Decrypt the message
        decrypted = rsa_decryption_with_d(encrypted, d, integer_N)
        print('Decrypted:', decrypted)
        print('Number of incorrect factorizations:', counter)
    else:
        # Increment the counter if the correct factors are not found
        counter += 1
except ValueError as e:
    print(e)