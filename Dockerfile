FROM python:3.12-bookworm

COPY --from=ghcr.io/astral-sh/uv:0.12 /uv /uvx /bin/

# venv は bind mount される /dr_tools の外に置く。
# container は host の UID で動き /opt/venv に書けないので、bytecode は build 時に作る。
ENV UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_PYTHON_DOWNLOADS=never \
    UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    PATH="/opt/venv/bin:${PATH}"

WORKDIR /dr_tools

COPY pyproject.toml uv.lock ./
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --extra dev --no-install-project

COPY . .
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --extra dev

CMD ["sleep", "infinity"]
