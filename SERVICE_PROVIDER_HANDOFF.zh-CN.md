# 给官网服务商的交付说明

我们希望把 DrugReflector 部署为官网 `www.tiniuni-bio.com` 的内置功能模块，而不是跳转到外部网站。

推荐访问路径：

```text
https://www.tiniuni-bio.com/drugreflector/
```

项目已经准备好官网内置部署包：

```text
deploy/official-site-module/
```

服务商只需要：

1. 把整个项目上传到服务器，例如 `/opt/drugreflector`
2. 确认 `checkpoints/` 下存在 3 个模型文件
3. 执行 Docker Compose 启动命令
4. 在官网 Nginx 里加入两个反向代理路径
5. 在官网后台新增“AI 药物筛选”入口，链接到 `/drugreflector/`

启动命令：

```bash
cd /opt/drugreflector
cp deploy/official-site-module/.env.example deploy/official-site-module/.env
docker compose --env-file deploy/official-site-module/.env -f deploy/official-site-module/docker-compose.yml up -d --build
```

官网 Nginx 配置文件参考：

```text
deploy/official-site-module/nginx/official-site-location.conf
```

部署说明全文：

```text
deploy/official-site-module/README.zh-CN.md
```

注意：

- 这是科研计算工具，不是普通静态页面，需要 Docker / Python / PyTorch 后端。
- 推荐服务器配置 4 核 CPU、16 GB 内存起步，磁盘 50 GB 以上。
- 必须使用 HTTPS。
- 官网需要允许至少 300 MB 上传。
