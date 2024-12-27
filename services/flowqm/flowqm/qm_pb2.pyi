from typing import ClassVar as _ClassVar
from typing import Iterable as _Iterable
from typing import Mapping as _Mapping
from typing import Optional as _Optional
from typing import Union as _Union

from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from google.protobuf.internal import containers as _containers

DESCRIPTOR: _descriptor.FileDescriptor

class Cell(_message.Message):
    __slots__ = ["values"]
    VALUES_FIELD_NUMBER: _ClassVar[int]
    values: _containers.RepeatedScalarFieldContainer[float]
    def __init__(self, values: _Optional[_Iterable[float]] = ...) -> None: ...

class Atoms(_message.Message):
    __slots__ = ["elements", "coordinates", "charge", "multiplicity", "pbc", "cell"]
    ELEMENTS_FIELD_NUMBER: _ClassVar[int]
    COORDINATES_FIELD_NUMBER: _ClassVar[int]
    CHARGE_FIELD_NUMBER: _ClassVar[int]
    MULTIPLICITY_FIELD_NUMBER: _ClassVar[int]
    PBC_FIELD_NUMBER: _ClassVar[int]
    CELL_FIELD_NUMBER: _ClassVar[int]
    elements: _containers.RepeatedScalarFieldContainer[str]
    coordinates: _containers.RepeatedScalarFieldContainer[float]
    charge: int
    multiplicity: int
    pbc: bool
    cell: Cell
    def __init__(
        self,
        elements: _Optional[_Iterable[str]] = ...,
        coordinates: _Optional[_Iterable[float]] = ...,
        charge: _Optional[int] = ...,
        multiplicity: _Optional[int] = ...,
        pbc: bool = ...,
        cell: _Optional[_Union[Cell, _Mapping]] = ...,
    ) -> None: ...

class CalculationParams(_message.Message):
    __slots__ = ["calc_type", "theory_level", "basis_set", "additional_params"]
    class AdditionalParamsEntry(_message.Message):
        __slots__ = ["key", "value"]
        KEY_FIELD_NUMBER: _ClassVar[int]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        key: str
        value: str
        def __init__(
            self, key: _Optional[str] = ..., value: _Optional[str] = ...
        ) -> None: ...

    CALC_TYPE_FIELD_NUMBER: _ClassVar[int]
    THEORY_LEVEL_FIELD_NUMBER: _ClassVar[int]
    BASIS_SET_FIELD_NUMBER: _ClassVar[int]
    ADDITIONAL_PARAMS_FIELD_NUMBER: _ClassVar[int]
    calc_type: str
    theory_level: str
    basis_set: str
    additional_params: _containers.ScalarMap[str, str]
    def __init__(
        self,
        calc_type: _Optional[str] = ...,
        theory_level: _Optional[str] = ...,
        basis_set: _Optional[str] = ...,
        additional_params: _Optional[_Mapping[str, str]] = ...,
    ) -> None: ...

class QMRequest(_message.Message):
    __slots__ = ["atoms", "params"]
    ATOMS_FIELD_NUMBER: _ClassVar[int]
    PARAMS_FIELD_NUMBER: _ClassVar[int]
    atoms: Atoms
    params: CalculationParams
    def __init__(
        self,
        atoms: _Optional[_Union[Atoms, _Mapping]] = ...,
        params: _Optional[_Union[CalculationParams, _Mapping]] = ...,
    ) -> None: ...

class QMResponse(_message.Message):
    __slots__ = ["success", "error_message", "atoms", "energy", "additional_data"]
    SUCCESS_FIELD_NUMBER: _ClassVar[int]
    ERROR_MESSAGE_FIELD_NUMBER: _ClassVar[int]
    ATOMS_FIELD_NUMBER: _ClassVar[int]
    ENERGY_FIELD_NUMBER: _ClassVar[int]
    ADDITIONAL_DATA_FIELD_NUMBER: _ClassVar[int]
    success: bool
    error_message: str
    atoms: Atoms
    energy: float
    additional_data: str
    def __init__(
        self,
        success: bool = ...,
        error_message: _Optional[str] = ...,
        atoms: _Optional[_Union[Atoms, _Mapping]] = ...,
        energy: _Optional[float] = ...,
        additional_data: _Optional[str] = ...,
    ) -> None: ...
