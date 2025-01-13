import streamlit as st
import yaml
from kubernetes import client, config
from kubernetes.client.rest import ApiException


def load_kube_config():
    try:
        config.load_kube_config()
    except Exception:
        config.load_incluster_config()


def create_namespace_if_not_exists(namespace: str, logs=None):
    if logs is None:
        logs = []

    core_v1 = client.CoreV1Api()
    try:
        core_v1.read_namespace(name=namespace)
        logs.append(f"Namespace '{namespace}' already exists.")
    except ApiException as e:
        if e.status == 404:
            logs.append(f"Namespace '{namespace}' not found. Creating it...")
            ns_body = client.V1Namespace(metadata=client.V1ObjectMeta(name=namespace))
            try:
                core_v1.create_namespace(ns_body)
                logs.append(f"Namespace '{namespace}' created.")
            except ApiException as e2:
                logs.append(f"[ERROR] Failed to create Namespace '{namespace}': {e2}")
        else:
            logs.append(f"[ERROR] {e}")

    return logs


def get_flowapp_url(namespace: str, logs=None) -> str:
    if logs is None:
        logs = []

    core_v1 = client.CoreV1Api()

    try:
        svc = core_v1.read_namespaced_service("flowapp", namespace)
    except ApiException as e:
        logs.append(f"[ERROR] Could not read Service 'flowapp': {e}")
        return ""

    if not svc.spec.ports or not svc.spec.ports[0].node_port:
        logs.append("[ERROR] flowapp Service is not NodePort or port not found.")
        return ""

    node_port = svc.spec.ports[0].node_port

    node_list = core_v1.list_node()
    if not node_list.items:
        logs.append("[ERROR] No nodes found in the cluster.")
        return ""

    node_ip = None
    for node in node_list.items:
        addresses = node.status.addresses or []
        external_ip = next(
            (a.address for a in addresses if a.type == "ExternalIP"), None
        )
        internal_ip = next(
            (a.address for a in addresses if a.type == "InternalIP"), None
        )

        if external_ip:
            node_ip = external_ip
            break
        elif internal_ip:
            node_ip = internal_ip
            break

    if not node_ip:
        logs.append("[ERROR] Could not find any ExternalIP or InternalIP in nodes.")
        return ""

    url = f"http://{node_ip}:{node_port}"
    return url


def create_resources_from_file(
    file_path, namespace="default", create_flowqm=False, logs=None
):
    if logs is None:
        logs = []

    with open(file_path, "r") as f:
        manifests = list(yaml.safe_load_all(f))

    apps_v1 = client.AppsV1Api()
    core_v1 = client.CoreV1Api()

    for manifest in manifests:
        kind = manifest.get("kind")
        name = manifest["metadata"]["name"]

        if name == "flowqm" and not create_flowqm:
            logs.append(f"Skipping {kind} '{name}' because flowqm is not selected.")
            continue
        elif name == "flowapp":
            pass
        else:
            logs.append(f"Skipping unknown {kind} '{name}'")
            continue

        if kind == "Deployment":
            logs.append(f"Creating Deployment '{name}' in '{namespace}' ...")
            try:
                resp = apps_v1.create_namespaced_deployment(
                    namespace=namespace, body=manifest
                )
                logs.append(
                    f"Deployment '{name}' created. (replicas={resp.spec.replicas})"
                )
            except ApiException as e:
                if e.status == 409:
                    logs.append(
                        f"Deployment '{name}' already exists. Skipping creation."
                    )
                else:
                    logs.append(f"[ERROR] {e}")

        elif kind == "Service":
            logs.append(f"Creating Service '{name}' in '{namespace}' ...")
            try:
                core_v1.create_namespaced_service(namespace=namespace, body=manifest)
                logs.append(f"Service '{name}' created.")
            except ApiException as e:
                if e.status == 409:
                    logs.append(f"Service '{name}' already exists. Skipping creation.")
                else:
                    logs.append(f"[ERROR] {e}")
        else:
            logs.append(f"Skipping unknown kind '{kind}'")

    logs.append("Finished creating selected resources.")
    return logs


def delete_resources_from_file(
    file_path, namespace="default", delete_flowqm=False, logs=None
):
    if logs is None:
        logs = []

    with open(file_path, "r") as f:
        manifests = list(yaml.safe_load_all(f))

    apps_v1 = client.AppsV1Api()
    core_v1 = client.CoreV1Api()

    for manifest in manifests:
        if manifest.get("kind") == "Service":
            name = manifest["metadata"]["name"]
            if name == "flowqm" and not delete_flowqm:
                logs.append(f"Skipping Service '{name}' (flowqm not selected).")
                continue
            elif name == "flowapp":
                pass
            else:
                logs.append(f"Skipping unknown Service '{name}'")
                continue

            logs.append(f"Deleting Service '{name}' from '{namespace}' ...")
            try:
                core_v1.delete_namespaced_service(name=name, namespace=namespace)
                logs.append(f"Service '{name}' deleted.")
            except ApiException as e:
                if e.status == 404:
                    logs.append(f"Service '{name}' not found. Skipping.")
                else:
                    logs.append(f"[ERROR] {e}")

    for manifest in manifests:
        if manifest.get("kind") == "Deployment":
            name = manifest["metadata"]["name"]
            if name == "flowqm" and not delete_flowqm:
                logs.append(f"Skipping Deployment '{name}' (flowqm not selected).")
                continue
            elif name == "flowapp":
                pass
            else:
                logs.append(f"Skipping unknown Deployment '{name}'")
                continue

            logs.append(f"Deleting Deployment '{name}' from '{namespace}' ...")
            try:
                apps_v1.delete_namespaced_deployment(name=name, namespace=namespace)
                logs.append(f"Deployment '{name}' deleted.")
            except ApiException as e:
                if e.status == 404:
                    logs.append(f"Deployment '{name}' not found. Skipping.")
                else:
                    logs.append(f"[ERROR] {e}")

    logs.append("Finished deleting selected resources.")
    return logs


def main():
    manifest_path = "manifest.yaml"
    st.title("FlowApp Apply/Delete")

    # ユーザ入力欄
    namespace = st.text_input("Namespace to use:", value="default")
    flowapp_checkbox = st.checkbox("flowapp (required)", value=True, disabled=True)  # noqa: F841
    flowqm_checkbox = st.checkbox(
        "flowqm (optional, QM calculation of psi4)", value=False
    )
    flowapp_url = None
    logs = []

    col1, col2 = st.columns(2)

    with col1:
        if st.button("Create", use_container_width=True, type="primary"):
            load_kube_config()
            logs.append("[ACTION] Create resources:")

            logs = create_namespace_if_not_exists(namespace, logs=logs)

            try:
                logs = create_resources_from_file(
                    file_path=manifest_path,
                    namespace=namespace,
                    create_flowqm=flowqm_checkbox,
                    logs=logs,
                )
            except Exception as e:
                logs.append(f"[ERROR] {e}")

            flowapp_url = get_flowapp_url(namespace, logs=logs)

    with col2:
        if st.button("Delete", use_container_width=True):
            load_kube_config()
            logs.append("[ACTION] Delete resources:")

            try:
                logs = delete_resources_from_file(
                    file_path=manifest_path,
                    namespace=namespace,
                    delete_flowqm=flowqm_checkbox,
                    logs=logs,
                )
            except Exception as e:
                logs.append(f"[ERROR] {e}")

    st.text_area("Logs:", value="\n".join(logs), height=300)
    if flowapp_url:
        logs.append(f"flowapp URL: {flowapp_url}")
        st.link_button(
            "Access flowapp",
            flowapp_url,
            type="primary",
            use_container_width=True,
        )


if __name__ == "__main__":
    main()
