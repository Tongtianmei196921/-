# DrugReflector 完整线上部署方案

当前项目分成两层：

- `Vercel`：只部署 React 前端网页。
- `Docker 后端`：部署 FastAPI + PyTorch + DrugReflector checkpoint，用于真实预测、GEO 解析、报告生成。

已经部署好的前端地址：

```text
https://drugreflector.vercel.app
```

注意：Vercel 上的前端不能单独完成预测，因为模型推理需要 Python、PyTorch、Scanpy 和 3 个模型 checkpoint。后端需要部署到支持 Docker 和较大内存的云平台。

## 推荐方案

优先推荐：

1. `Railway` 或 `Render`
2. `Hugging Face Spaces Docker`
3. 普通 Linux 云服务器 + Docker

`Railway/Render` 更像正式 Web API；`Hugging Face Spaces` 更像科研 demo；Linux 云服务器最可控。

## 后端部署要求

后端至少需要：

- Docker 构建支持
- 2 CPU 或以上
- 4GB 内存起步，8GB 更稳
- 能访问外网下载 Zenodo checkpoint
- 公网 HTTPS URL

当前 `Dockerfile` 会在构建时自动运行：

```bash
zenodo_get --output-dir checkpoints 16912444
```

所以不需要把 `.pt` 模型文件提交到 GitHub。

## Railway / Render 部署思路

1. 在平台中新建 Docker Web Service。
2. 连接 GitHub 仓库：

```text
https://github.com/Tongtianmei196921/-
```

3. 选择从根目录 Dockerfile 构建。
4. 暴露端口使用平台自动注入的 `PORT`。
5. 部署成功后确认：

```text
https://你的后端域名/health
https://你的后端域名/api/checkpoints
```

`/api/checkpoints` 应该显示 3 个模型都是 `true`。

## 把 Vercel 前端接到后端

后端部署成功后，在 Vercel 项目里设置环境变量：

```text
VITE_API_BASE_URL=https://你的后端域名
```

然后重新部署 Vercel。

重新部署后，前端会把这些请求发到公网后端：

- `/api/checkpoints`
- `/api/predict`
- `/api/geo/preview`
- `/api/geo/predict`
- `/api/prepare-input`
- `/api/report/docx`

## 本地验证

后端本地运行：

```bash
docker build -t drugreflector-api .
docker run --rm -p 8000:8000 drugreflector-api
```

前端本地连接远程后端：

```bash
cd frontend
set VITE_API_BASE_URL=https://你的后端域名
npm run build
```

## 重要边界

- Vercel 适合前端，不适合直接跑本项目完整后端。
- 线上预测会消耗内存和 CPU，GEO 数据越大越慢。
- 免费平台可能会休眠，首次访问可能较慢。
- 如果给客户正式使用，建议使用付费 Docker Web Service 或独立云服务器。
