#  baxter
Bacterial variant exercise generator

## Building and starting docker services


```bash
export MYSQL_ROOT_PASSWORD=something
docker compose up -d baxter-mariadb
```
Docker-compose reads `docker-compose.yml` file.
Note that `-d` will put the service  in the background.


```bash
export MYSQL_ROOT_PASSWORD=something
docker compose up -d baxter-backend
```

```bash
export MYSQL_ROOT_PASSWORD=something
docker compose up -d baxter-frontend
```

Note that the sub-commands `create` and `start` can be sued to decouple these two operations:
```bash
docker compose create baxter-frontend
docker compose start baxter-backend
```

## Troubleshooting

```bash
bash> docker compose up -d baxter-frontend
[+] Running 0/1
 ⠿ baxter-frontend Error                                                                                                                       1.7s
[+] Building 0.0s (0/0)                                                                                                                             
failed to solve: rpc error: code = Unavailable desc = connection error: 
desc = "transport: Error while dialing unable to upgrade to h2c, received 404"
```
is a meaningless error message. It could be anything. Try running dovker build in the
corresponding directory to see the output and th actual problem
