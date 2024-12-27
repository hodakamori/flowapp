import threading
import time
from concurrent import futures

import grpc
import pytest
from flowqm import qm_pb2, qm_pb2_grpc
from flowqm.service import QMService


@pytest.fixture(scope="module")
def grpc_server():
    server = grpc.server(futures.ThreadPoolExecutor(max_workers=4))
    qm_pb2_grpc.add_QMServiceServicer_to_server(QMService(), server)
    port = 50051
    server.add_insecure_port(f"[::]:{port}")

    server_thread = threading.Thread(target=server.start)
    server_thread.start()

    time.sleep(1)

    yield f"localhost:{port}"

    server.stop(0)
    server_thread.join()


def test_qm_service(grpc_server) -> None:
    with grpc.insecure_channel("localhost:50051") as channel:
        stub = qm_pb2_grpc.QMServiceStub(channel)
        request = qm_pb2.QMRequest(
            atoms=qm_pb2.Atoms(
                elements=["H", "O", "H"],
                coordinates=[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                charge=0,
                multiplicity=1,
            ),
            params=qm_pb2.CalculationParams(
                calc_type="opt",
                theory_level="hf",
                basis_set="sto-3g",
                additional_params={"key1": "value1"},
            ),
        )
        response = stub.CalculateQM(request)
        assert response.success
        assert response.error_message == ""
        assert response.atoms.elements == ["H", "O", "H"]
        assert len(response.atoms.coordinates) == 3
