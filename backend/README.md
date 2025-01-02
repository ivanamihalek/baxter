

```bash
python3 -m venv venv
. venv/bin/activate
pip install -r requirements.txt
```



```bash
docker  build -t baxter-backend-img --load  --build-arg mysql_password=fhghgqieu .
docker run -d --name baxter-backend-ctnr -p 8000:8000 baxter-backend-img
```

