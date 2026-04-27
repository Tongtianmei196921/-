# 自动更新说明

当前线上系统分为两部分：

- 前端：Vercel，地址 `https://drugreflector.vercel.app`
- 后端：Hugging Face Space，地址 `https://qiangaoqing-drugreflector-api.hf.space`

仓库已经配置 GitHub Actions。以后本地修改代码后，只需要：

```bash
git add .
git commit -m "更新说明"
git push
```

GitHub 会自动根据改动内容部署对应平台。

## 需要在 GitHub 里配置一次的 Secrets

进入仓库：

```text
https://github.com/Tongtianmei196921/-
```

打开：

```text
Settings -> Secrets and variables -> Actions -> New repository secret
```

添加两个 secret：

```text
VERCEL_TOKEN
HF_TOKEN
```

## VERCEL_TOKEN 获取方式

1. 打开 Vercel。
2. 进入 Account Settings。
3. 找到 Tokens。
4. 创建一个新 token。
5. 复制到 GitHub Secret：`VERCEL_TOKEN`。

项目 ID 和组织 ID 已经写在工作流里：

```text
VERCEL_ORG_ID=team_qlUVYwRm9yI7cSs7CrearCAQ
VERCEL_PROJECT_ID=prj_wVmJTNXhu85s1LEsgudazZ9M5T6l
```

## HF_TOKEN 获取方式

1. 打开 Hugging Face。
2. 进入 Settings -> Access Tokens。
3. 创建一个有写入权限的 token。
4. 复制到 GitHub Secret：`HF_TOKEN`。

后端 Space 已经固定为：

```text
Qiangaoqing/drugreflector-api
```

## 什么时候会自动部署前端

改动这些文件并推送到 `main` 后，会自动部署 Vercel：

- `frontend/**`
- `vercel.json`
- `.github/workflows/deploy-frontend-vercel.yml`

## 什么时候会自动部署后端

改动这些文件并推送到 `main` 后，会自动部署 Hugging Face Space：

- `Dockerfile`
- `pyproject.toml`
- `hf-space.README.md`
- `drugreflector/**`
- `signature_refinement/**`
- `.github/workflows/deploy-backend-huggingface.yml`

## 自动更新后的验证

前端：

```text
https://drugreflector.vercel.app
```

后端健康检查：

```text
https://qiangaoqing-drugreflector-api.hf.space/health
```

后端模型检查：

```text
https://qiangaoqing-drugreflector-api.hf.space/api/checkpoints
```

`all_ready` 应为 `true`。

## 注意

- GitHub Actions 只有配置好 `VERCEL_TOKEN` 和 `HF_TOKEN` 后才能真正自动部署。
- 免费 Hugging Face Space 可能会休眠，首次访问会较慢。
- 如果只改 README 或本地文档，不一定会触发部署。
