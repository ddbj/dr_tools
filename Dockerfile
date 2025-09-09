FROM python:3.12-bookworm

# パッケージリストの更新とPythonのインストール
RUN apt update && \
    apt install -y \
    curl \
    less && \
    rm -rf /var/lib/apt/lists/* && \
    ln -s /usr/bin/python3 /usr/bin/python

RUN alias ll='ls -l'

# biopythonのインストール
RUN python3 -m pip install --break-system-packages biopython==1.83

WORKDIR /dr_tools
COPY . .
RUN python3 -m pip install -e .[dev]

ENTRYPOINT [ "" ]
# Pythonのバージョン確認用コマンド
CMD ["python3", "--version"]
