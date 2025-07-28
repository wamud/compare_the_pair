# With modifications by Ant

# Copyright 2022 Oscar Higgott

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#      http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import stim
from typing import Callable, Set, List, Dict, Tuple, Optional
from dataclasses import dataclass
import math


def append_anti_basis_error(circuit: stim.Circuit, targets: List[int], p: float, basis: str) -> None:
    if p > 0:
        if basis == "X":
            circuit.append_operation("Z_ERROR", targets, p)
        else:
            circuit.append_operation("X_ERROR", targets, p)


@dataclass
class CircuitGenParameters:
    code_name: str
    task: str
    rounds: int
    distance: int = None
    x_distance: int = None
    z_distance: int = None
    after_clifford_depolarization: float = 0
    before_round_data_depolarization: float = 0
    before_measure_flip_probability: float = 0
    after_reset_flip_probability: float = 0
    exclude_other_basis_detectors: bool = False

    def append_begin_round_tick(
            self,
            circuit: stim.Circuit,
            data_qubits: List[int]
    ) -> None:
        circuit.append_operation("TICK", [])
        if self.before_round_data_depolarization > 0:
            circuit.append_operation("DEPOLARIZE1", data_qubits, self.before_round_data_depolarization)
            
    def append_unitary_1(
            self,
            circuit: stim.Circuit,
            name: str,
            targets: List[int]
    ) -> None:
        circuit.append_operation(name, targets)
        
        if self.after_clifford_depolarization > 0:
            circuit.append_operation("DEPOLARIZE1", targets, self.after_clifford_depolarization) # comment out to take off error after Hadamard 
            pass
        
    def append_unitary_2(
            self,
            circuit: stim.Circuit,
            name: str,
            targets: List[int],
    ) -> None:
        circuit.append_operation(name, targets)
        if self.after_clifford_depolarization > 0:
            circuit.append_operation("DEPOLARIZE2", targets, self.after_clifford_depolarization)

    def append_reset(
            self,
            circuit: stim.Circuit,
            targets: List[int],
            basis: str = "Z"
    ) -> None:
        circuit.append_operation("R" + basis, targets)
        append_anti_basis_error(circuit, targets, self.after_reset_flip_probability, basis)

    def append_measure(self, circuit: stim.Circuit, targets: List[int], basis: str = "Z") -> None:
        append_anti_basis_error(circuit, targets, self.before_measure_flip_probability, basis)
        circuit.append_operation("M" + basis, targets)

    def append_measure_reset(
            self,
            circuit: stim.Circuit,
            targets: List[int],
            basis: str = "Z"
    ) -> None:
        append_anti_basis_error(circuit, targets, self.before_measure_flip_probability, basis)
        circuit.append_operation("MR" + basis, targets)
        append_anti_basis_error(circuit, targets, self.after_reset_flip_probability, basis)


def finish_surface_code_circuit(
        coord_to_index: Callable[[complex], int],
        data_coords: Set[complex],
        x_measure_coords: Set[complex],
        z_measure_coords: Set[complex],
        params: CircuitGenParameters,
        x_order: List[complex],
        z_order: List[complex],
        x_observable: List[complex],
        z_observable: List[complex],
        is_memory_x: bool,
        *,
        exclude_other_basis_detectors: bool = False,
        wraparound_length: Optional[int] = None
) -> stim.Circuit:
    if params.rounds < 1:
        raise ValueError("Need rounds >= 1")
    if params.distance is not None and params.distance < 2:
        raise ValueError("Need a distance >= 2")
    if params.x_distance is not None and (params.x_distance < 2 or
                                          params.z_distance < 2):
        raise ValueError("Need a distance >= 2")

    chosen_basis_observable = x_observable if is_memory_x else z_observable
    chosen_basis_measure_coords = x_measure_coords if is_memory_x else z_measure_coords

    # Index the measurement qubits and data qubits.
    p2q: Dict[complex, int] = {}
    for q in data_coords:
        p2q[q] = coord_to_index(q)

    for q in x_measure_coords:
        p2q[q] = coord_to_index(q)

    for q in z_measure_coords:
        p2q[q] = coord_to_index(q)

    q2p: Dict[int, complex] = {v: k for k, v in p2q.items()}

    data_qubits = [p2q[q] for q in data_coords]
    measurement_qubits = [p2q[q] for q in x_measure_coords]
    measurement_qubits += [p2q[q] for q in z_measure_coords]
    x_measurement_qubits = [p2q[q] for q in x_measure_coords]

    all_qubits: List[int] = []
    all_qubits += data_qubits + measurement_qubits

    all_qubits.sort()
    data_qubits.sort()
    measurement_qubits.sort()
    x_measurement_qubits.sort()

    # Reverse index the measurement order used for defining detectors
    data_coord_to_order: Dict[complex, int] = {}
    measure_coord_to_order: Dict[complex, int] = {}
    for q in data_qubits:
        data_coord_to_order[q2p[q]] = len(data_coord_to_order)
    for q in measurement_qubits:
        measure_coord_to_order[q2p[q]] = len(measure_coord_to_order)

    # List out CNOT gate targets using given interaction orders.
    cnot_targets: List[List[int]] = [[], [], [], []]
    for k in range(4):
        for measure in sorted(x_measure_coords, key=lambda c: (c.real, c.imag)):
            data = measure + x_order[k]
            if data in p2q:
                cnot_targets[k].append(p2q[measure])
                cnot_targets[k].append(p2q[data])
            elif wraparound_length is not None:
                data_wrapped = (data.real % wraparound_length) + (data.imag % wraparound_length) * 1j
                cnot_targets[k].append(p2q[measure])
                cnot_targets[k].append(p2q[data_wrapped])

        for measure in sorted(z_measure_coords, key=lambda c: (c.real, c.imag)):
            data = measure + z_order[k]
            if data in p2q:
                cnot_targets[k].append(p2q[data])
                cnot_targets[k].append(p2q[measure])
            elif wraparound_length is not None:
                data_wrapped = (data.real % wraparound_length) + (data.imag % wraparound_length) * 1j
                cnot_targets[k].append(p2q[data_wrapped])
                cnot_targets[k].append(p2q[measure])

    # Build the repeated actions that make up the surface code cycle
    cycle_actions = stim.Circuit()
    params.append_begin_round_tick(cycle_actions, data_qubits)
    params.append_unitary_1(cycle_actions, "H", x_measurement_qubits)
    for targets in cnot_targets:
        cycle_actions.append_operation("TICK", [])
        params.append_unitary_2(cycle_actions, "CNOT", targets)
    cycle_actions.append_operation("TICK", [])
    params.append_unitary_1(cycle_actions, "H", x_measurement_qubits)
    cycle_actions.append_operation("TICK", [])
    params.append_measure_reset(cycle_actions, measurement_qubits)

    # Build the start of the circuit, getting a state that's ready to cycle
    # In particular, the first cycle has different detectors and so has to be handled special.
    head = stim.Circuit()
    for k, v in sorted(q2p.items()):
        head.append_operation("QUBIT_COORDS", [k], [v.real, v.imag])
    params.append_reset(head, data_qubits, "ZX"[is_memory_x])
    params.append_reset(head, measurement_qubits)
    head += cycle_actions
    for measure in sorted(chosen_basis_measure_coords, key=lambda c: (c.real, c.imag)):
        head.append_operation(
            "DETECTOR",
            [stim.target_rec(-len(measurement_qubits) + measure_coord_to_order[measure])],
            [measure.real, measure.imag, 0.0]
        )

    # Build the repeated body of the circuit, including the detectors comparing to previous cycles.
    body = cycle_actions.copy()
    m = len(measurement_qubits)
    body.append_operation("SHIFT_COORDS", [], [0.0, 0.0, 1.0])
    for m_index in measurement_qubits:
        m_coord = q2p[m_index]
        k = len(measurement_qubits) - measure_coord_to_order[m_coord] - 1
        if not exclude_other_basis_detectors or m_coord in chosen_basis_measure_coords:
            body.append_operation(
                "DETECTOR",
                [stim.target_rec(-k - 1), stim.target_rec(-k - 1 - m)],
                [m_coord.real, m_coord.imag, 0.0]
            )

    # Build the end of the circuit, getting out of the cycle state and terminating.
    # In particular, the data measurements create detectors that have to be handled special.
    # Also, the tail is responsible for identifying the logical observable.
    tail = stim.Circuit()
    params.append_measure(tail, data_qubits, "ZX"[is_memory_x])
    # Detectors
    for measure in sorted(chosen_basis_measure_coords, key=lambda c: (c.real, c.imag)):
        detectors: List[int] = []
        for delta in z_order:
            data = measure + delta
            if data in p2q:
                detectors.append(-len(data_qubits) + data_coord_to_order[data])
            elif wraparound_length is not None:
                data_wrapped = (data.real % wraparound_length) + (data.imag % wraparound_length) * 1j
                detectors.append(-len(data_qubits) + data_coord_to_order[data_wrapped])
        detectors.append(-len(data_qubits) - len(measurement_qubits) + measure_coord_to_order[measure])
        detectors.sort(reverse=True)
        tail.append_operation("DETECTOR", [stim.target_rec(x) for x in detectors], [measure.real, measure.imag, 1.0])

    # Logical observable
    obs_inc: List[int] = []
    for q in chosen_basis_observable:
        obs_inc.append(-len(data_qubits) + data_coord_to_order[q])
    obs_inc.sort(reverse=True)
    tail.append_operation("OBSERVABLE_INCLUDE", [stim.target_rec(x) for x in obs_inc], 0.0)

    # Combine to form final circuit.
    return head + body * (params.rounds - 1) + tail


def generate_rotated_surface_code_circuit(
        params: CircuitGenParameters,
        is_memory_x: bool,
        x_cnot_order,
        z_cnot_order
) -> stim.Circuit:
    if params.distance is not None:
        x_distance = params.distance
        z_distance = params.distance
    else:
        x_distance = params.x_distance
        z_distance = params.z_distance

    # Place data qubits
    data_coords: Set[complex] = set()
    x_observable: List[complex] = []
    z_observable: List[complex] = []
    for x in [i + 0.5 for i in range(z_distance)]:
        for y in [i + 0.5 for i in range(x_distance)]:
            q = x * 2 + y * 2 * 1j
            data_coords.add(q)
            if y == 0.5:
                z_observable.append(q)
            if x == 0.5:
                x_observable.append(q)

    # Place measurement qubits.
    x_measure_coords: Set[complex] = set()
    z_measure_coords: Set[complex] = set()
    for x in range(z_distance + 1):
        for y in range(x_distance + 1):
            q = x * 2 + y * 2j
            on_boundary_1 = x == 0 or x == z_distance
            on_boundary_2 = y == 0 or y == x_distance
            parity = (x % 2) != (y % 2)
            if on_boundary_1 and parity:
                continue
            if on_boundary_2 and not parity:
                continue
            if parity:
                x_measure_coords.add(q)
            else:
                z_measure_coords.add(q)

    # Default interaction orders so that hook errors run against the error grain instead of with it:
    x_order: List[complex] = [1 + 1j, -1 + 1j, 1 - 1j, -1 - 1j]   # d: this order is 2310 (counting 0 to 3 clockwise from top left)
    z_order: List[complex] = [1 + 1j, 1 - 1j, -1 + 1j, -1 - 1j]   # d: this order is 2130 

    # Ant: above is original ordering. An order of 0,1,2,3 (counting from top left of each stabiliser square in clockwise direction) is:
    clockwise_order: List[complex] = [-1 - 1j, 1 - 1j, 1 + 1j, -1 + 1j]
    # (This is the positions of a square's corners on an imaginary and real axis but with up being the negative imaginary direction)
    # Define new orders:
    x_order = [clockwise_order[i] for i in x_cnot_order]   # x hook error if last two CNOT's are 12, 21, 03 or 30 
    z_order = [clockwise_order[i] for i in z_cnot_order]   # z hook error if last two CNOT's are 01, 10, 23 or 32 

        
    def coord_to_idx(q: complex) -> int:
        q = q - math.fmod(q.real, 2) * 1j
        return int(q.real + q.imag * (z_distance + 0.5))

    return finish_surface_code_circuit(
        coord_to_idx,
        data_coords,
        x_measure_coords,
        z_measure_coords,
        params,
        x_order,
        z_order,
        x_observable,
        z_observable,
        is_memory_x,
        exclude_other_basis_detectors=params.exclude_other_basis_detectors
    )


def _generate_unrotated_surface_or_toric_code_circuit(
        params: CircuitGenParameters,
        is_memory_x: bool,
        is_toric: bool,
        z_cnot_order,
        x_cnot_order
) -> stim.Circuit:
    d = params.distance
    assert params.rounds > 0

    # Place qubits
    data_coords: Set[complex] = set()
    x_measure_coords: Set[complex] = set()
    z_measure_coords: Set[complex] = set()
    x_observable: List[complex] = []
    z_observable: List[complex] = []
    length = 2 * d if is_toric else 2 * d - 1
    for x in range(length):
        for y in range(length):
            q = x + y * 1j
            parity = (x % 2) != (y % 2)
            if parity:
                if x % 2 == 0:
                    z_measure_coords.add(q)
                else:
                    x_measure_coords.add(q)
            else:
                data_coords.add(q)
                if x == 0:
                    x_observable.append(q)
                if y == 0:
                    z_observable.append(q)

    # Define interaction order. Doesn't matter so much for unrotated.
    
    # order: List[complex] = [1, 1j, -1j, -1]  # Ant: this is the original order: 1203


    # Ant: [ above is original ordering. An order of 0,1,2,3 (counting from top of each square in clockwise direction) is:
    clockwise_order: List[complex] = [-1j, 1, 1j, -1]

    # Set CNOT orderings for unrotated code (previously "unrot_cnot_order")
    z_order = [clockwise_order[i] for i in z_cnot_order]
    x_order = [clockwise_order[i] for i in x_cnot_order]

    # ]

    def coord_to_idx(q: complex) -> int:
        return int(q.real + q.imag * length)

    # Delegate.
    return finish_surface_code_circuit(
        coord_to_idx,
        data_coords,
        x_measure_coords,
        z_measure_coords,
        params,
        x_order,
        z_order,
        x_observable,
        z_observable,
        is_memory_x,
        exclude_other_basis_detectors=params.exclude_other_basis_detectors,
        wraparound_length=2 * d if is_toric else None
    )


def generate_surface_or_toric_code_circuit_from_params(params: CircuitGenParameters,z_cnot_order,x_cnot_order) -> stim.Circuit:
    if params.code_name == "surface_code":
        if params.task == "rotated_memory_x":
            return generate_rotated_surface_code_circuit(params, True,z_cnot_order,x_cnot_order)
        elif params.task == "rotated_memory_z":
            return generate_rotated_surface_code_circuit(params, False,z_cnot_order,x_cnot_order)
        elif params.task == "unrotated_memory_x":
            if params.distance is None:
                raise NotImplementedError('Rectangular unrotated memories are '
                                          'not currently supported')
            return _generate_unrotated_surface_or_toric_code_circuit(
                params=params,
                is_memory_x=True,
                is_toric=False,
                z_cnot_order=z_cnot_order,
                x_cnot_order = x_cnot_order)
        elif params.task == "unrotated_memory_z":
            if params.distance is None:
                raise NotImplementedError('Rectangular unrotated memories are '
                                          'not currently supported')
            return _generate_unrotated_surface_or_toric_code_circuit(
                params=params,
                is_memory_x=False,
                is_toric=False,
                x_cnot_order=x_cnot_order,
                z_cnot_order=z_cnot_order)
    elif params.code_name == "toric_code":
        if params.distance is None:
            raise NotImplementedError('Rectangular toric codes are '
                                      'not currently supported')
        if params.task == "unrotated_memory_x":
            return _generate_unrotated_surface_or_toric_code_circuit(
                params=params,
                is_memory_x=True,
                is_toric=True,
                x_cnot_order=x_cnot_order,
                z_cnot_order=z_cnot_order)
        elif params.task == "unrotated_memory_z":
            return _generate_unrotated_surface_or_toric_code_circuit(
                params=params,
                is_memory_x=False,
                is_toric=True,
                x_cnot_order=x_cnot_order,
                z_cnot_order=z_cnot_order)

    raise ValueError(f"Unrecognised task: {params.task}")


def generate_circuit(
        code_task: str,
        *,
        rounds: int,
        distance: int = None,
        x_distance: int = None,
        z_distance: int = None,
        after_clifford_depolarization: float = 0.0,
        before_round_data_depolarization: float = 0.0,
        before_measure_flip_probability: float = 0.0,
        after_reset_flip_probability: float = 0.0,
        exclude_other_basis_detectors: bool = False,
        x_cnot_order: list = [2,3,1,0],
        z_cnot_order: list = [2,1,3,0]
) -> stim.Circuit:
    """Generates common circuits.

        The generated circuits can include configurable noise.

        The generated circuits include DETECTOR and OBSERVABLE_INCLUDE annotations so
        that their detection events and logical observables can be sampled.

        The generated circuits include TICK annotations to mark the progression of time.
        (E.g. so that converting them using `stimcirq.stim_circuit_to_cirq_circuit` will
        produce a `cirq.Circuit` with the intended moment structure.)

        Note that the toric_code circuits currently only include one of the two logical observables.

        Args:
            code_task: A string identifying the type of circuit to generate. Available
                code tasks are:
                    - "surface_code:rotated_memory_x"
                    - "surface_code:rotated_memory_z"
                    - "surface_code:unrotated_memory_x"
                    - "surface_code:unrotated_memory_z"
                    - "toric_code:unrotated_memory_x"
                    - "toric_code:unrotated_memory_z"
            distance: Defaults to None. The desired code distance of the generated
                circuit. The code distance is the minimum number of physical
                errors needed to cause a logical error. This parameter indirectly determines
                how many qubits the generated circuit uses.
            x_distance: Defaults to None. The desired code distance of the X
            logical
                operator in the generated circuit: the minimum number of X physical
                errors needed to cause a X logical error.
            z_distance: Defaults to None. The desired code distance of the Z
            logical
                operator in the generated circuit: the minimum number of Z physical
                errors needed to cause a Z logical error.
            rounds: How many times the measurement qubits in the generated circuit will
                be measured. Indirectly determines the duration of the generated
                circuit.
            after_clifford_depolarization: Defaults to 0. The probability (p) of
                `DEPOLARIZE1(p)` operations to add after every single-qubit Clifford
                operation and `DEPOLARIZE2(p)` operations to add after every two-qubit
                Clifford operation. The after-Clifford depolarizing operations are only
                included if this probability is not 0.
            before_round_data_depolarization: Defaults to 0. The probability (p) of
                `DEPOLARIZE1(p)` operations to apply to every data qubit at the start of
                a round of stabilizer measurements. The start-of-round depolarizing
                operations are only included if this probability is not 0.
            before_measure_flip_probability: Defaults to 0. The probability (p) of
                `X_ERROR(p)` operations applied to qubits before each measurement (X
                basis measurements use `Z_ERROR(p)` instead). The before-measurement
                flips are only included if this probability is not 0.
            after_reset_flip_probability: Defaults to 0. The probability (p) of
                `X_ERROR(p)` operations applied to qubits after each reset (X basis
                resets use `Z_ERROR(p)` instead). The after-reset flips are only
                included if this probability is not 0.
            exclude_other_basis_detectors: Defaults to False. If True, do not add
                detectors to measurement qubits that are measured in the opposite
                basis to the chosen basis of the logical observable.
            x_cnot_order: Defaults to [2,3,1,0] (same as stim's default for rotated code.
                For unrotated the default is 1203). This defines the order of the 
                CNOT gates for each X stabiliser measurement in the surface
                code. When looking at the diagram from 
                <generated_circuit>.diagram("detslice-with-ops-svg") each auxiliary 
                qubit is at the centre of a square of data qubits. The order counts
                0 to 3 in a clockwise direction from the top left data qubit in the rotated
                code, the top data qubit in the unrotated.
            z_cnot_order: Defaults to [2,1,3,0] (same as stim's default for rotated code.
                For unrotated the default is the same as unrotated x_cnot order, 1203). 
                This defines the order of the CNOT gates for each Z stabiliser 
                measurement in the surface code. When looking at the diagram from 
                <generated_circuit>.diagram("detslice-with-ops-svg") each auxiliary 
                qubit is at the centre of a square of data qubits. The order counts
                0 to 3 in a clockwise direction from the top left data qubit in the rotated
                code, the top data qubit in the unrotated.

            # I believe generating circuits with this code results in a different strong_ID to when they were generated by oscAntgencricuitv1.py even if the CNOT orders are the same. 


        Returns:
            The generated circuit.
        """
    if distance is not None:
        pass
    elif x_distance is not None and z_distance is not None:
        pass
    else:
        raise ValueError('Either the distance parameter or x_distance and '
                         'z_distance parameters must be specified')
    code_name, task = code_task.split(":")
    if code_name in ["surface_code", "toric_code"]:
        params = CircuitGenParameters(
            code_name=code_name,
            task=task,
            rounds=rounds,
            distance=distance,
            x_distance=x_distance,
            z_distance=z_distance,
            after_clifford_depolarization=after_clifford_depolarization,
            before_round_data_depolarization=before_round_data_depolarization,
            before_measure_flip_probability=before_measure_flip_probability,
            after_reset_flip_probability=after_reset_flip_probability,
            exclude_other_basis_detectors=exclude_other_basis_detectors,
        )
        return generate_surface_or_toric_code_circuit_from_params(params,x_cnot_order,z_cnot_order)
    else:
        raise ValueError(f"Code name {code_name} not recognised")