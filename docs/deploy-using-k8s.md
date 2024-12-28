### Kubernetes Deployment Guide
#### Prerequisites
- **microk8s**: Installed and configured.
- **Docker**: Installed for building images.
- **kubectl**: Configured to work with microk8s.

#### Steps

1. **Build Docker Images**
   ```bash
   docker build -t flowapp:latest .
   docker build -t flowqm:latest -f services/flowqm/Dockerfile .
   ```

2. **Import Images to microk8s**
   ```bash
   docker save flowapp:latest | microk8s.ctr image import -
   docker save flowqm:latest | microk8s.ctr image import -
   ```

3. **Deploy to Kubernetes**
   ```bash
   kubectl apply -f manifest.yaml
   ```

4. **Expose External Access**
   Edit the `flowapp` service to use `NodePort`:
   ```bash
   kubectl edit service flowapp
   ```
   Modify `spec.type` to `NodePort` and save changes.

5. **Access FlowApp**
   Find the NodePort:
   ```bash
   kubectl get services
   ```
   Access FlowApp in your browser:
   ```
   http://<HOST_IP>:<NODE_PORT>
   ```

---

## Cleanup
To remove the deployment:
```bash
microk8s.kubectl delete -f manifest.yaml
