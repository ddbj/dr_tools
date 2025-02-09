# Ubuntu 24.04 LTSをベースイメージとして使用
FROM ubuntu:24.04

# パッケージリストの更新とPythonのインストール
RUN apt-get update && \
    apt-get install -y \
    python3 \
    python3-pip \
    curl less \
    && rm -rf /var/lib/apt/lists/* && \
    ln -s /usr/bin/python3 /usr/bin/python

# biopythonのインストール
RUN pip3 install --break-system-packages biopython==1.83

RUN alias ll='ls -l'

# Pythonのバージョン確認用コマンド
CMD ["python3", "--version"] 