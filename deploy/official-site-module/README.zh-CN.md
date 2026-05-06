# DrugReflector 官网内置模块部署说明

这是给 `www.tiniuni-bio.com` 官网服务商使用的部署包。目标是把 DrugReflector 作为官网内置功能模块访问：

```text
https://www.tiniuni-bio.com/drugreflector/
```

服务商不需要理解模型算法，只需要按下面步骤部署。

## 一、服务器要求

- Linux 服务器，推荐 Ubuntu 22.04 / Debian 12 / CentOS 7+
- Docker 和 Docker Compose
- 建议 4 核 CPU、16 GB 内存起步
- 建议磁盘 50 GB 以上
- 官网必须支持 HTTPS
- 官网 Nginx / 网关需要允许 300 MB 上传

如果多人同时使用、大文件较多，建议 8 核 CPU、32 GB 内存。

## 二、一键启动

把整个项目目录上传到服务器，例如：

```text
/opt/drugreflector
```

进入目录：

```bash
cd /opt/drugreflector
cp deploy/official-site-module/.env.example deploy/official-site-module/.env
docker compose --env-file deploy/official-site-module/.env -f deploy/official-site-module/docker-compose.yml up -d --build
```

启动后先在服务器本机测试：

```bash
curl http://127.0.0.1:18080/drugreflector/
curl http://127.0.0.1:18080/drugreflector-api/api/checkpoints
```

第二个接口应返回 3 个模型文件均为 `true`。

## 三、模型文件

模型文件放在：

```text
checkpoints/
  model_fold_0.pt
  model_fold_1.pt
  model_fold_2.pt
```

如果交付包里已经带有这 3 个文件，服务商不需要处理。

如果没有模型文件，可以执行：

```bash
sh deploy/official-site-module/scripts/download-checkpoints.sh
```

Windows 测试环境可以执行：

```powershell
.\deploy\official-site-module\scripts\download-checkpoints.ps1
```

## 四、接入官网 Nginx

把下面配置加入 `www.tiniuni-bio.com` 的 HTTPS server block：

```nginx
client_max_body_size 300m;
proxy_read_timeout 600s;
proxy_send_timeout 600s;

location ^~ /drugreflector/ {
    proxy_pass http://127.0.0.1:18080/drugreflector/;
    proxy_http_version 1.1;
    proxy_set_header Host $host;
    proxy_set_header X-Real-IP $remote_addr;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
    proxy_set_header X-Forwarded-Proto $scheme;
}

location ^~ /drugreflector-api/ {
    proxy_pass http://127.0.0.1:18080/drugreflector-api/;
    proxy_http_version 1.1;
    proxy_set_header Host $host;
    proxy_set_header X-Real-IP $remote_addr;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
    proxy_set_header X-Forwarded-Proto $scheme;
}
```

配置完成后重载 Nginx：

```bash
nginx -t && systemctl reload nginx
```

## 五、官网入口

请在官网后台新增一个菜单或按钮：

```text
AI 药物筛选
```

链接：

```text
/drugreflector/
```

## 六、验收清单

- `https://www.tiniuni-bio.com/drugreflector/` 可以打开
- `https://www.tiniuni-bio.com/drugreflector-api/api/checkpoints` 返回 3 个模型均 ready
- 可以上传示例 CSV
- 可以运行预测
- 可以下载 SVG / PNG / CSV / Word 报告
- 手机端可以正常打开页面
- HTTPS 证书正常，不出现浏览器安全警告

## 七、数据和合规提醒

- 用户上传的数据只用于本次分析，不应公开展示
- 不建议把用户上传文件长期保存
- 后续如果作为正式对外服务，请让官网服务商确认 ICP / 公安备案覆盖范围
