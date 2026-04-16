#!/bin/bash

#@ Network.hppに関するヘッダーファイルは，このプリコンパイル`pch.hpp`に含めていない
#@ 変更があまりされないものだけを，`pch.hpp`に含めている

SCRIPT_DIR=$(cd $(dirname $0) && pwd)

# 必要な場合は以下を実行
(
   cd "${SCRIPT_DIR}/tetgen1.6.0"
   sh clean
   cmake -DCMAKE_BUILD_TYPE=Release ./
   make
)

# OSの判定
OS_TYPE=$(uname)

if [ "$OS_TYPE" = "Darwin" ]; then
   # macOSの場合
   echo "macOS detected"
   (
      cd "${SCRIPT_DIR}"
      /opt/homebrew/bin/g++-14 -v -std=gnu++23 -x c++-header -fopenmp -fconcepts -pthread -Ofast -march=native \
         -I"${SCRIPT_DIR}" \
         -I"${SCRIPT_DIR}/../include" \
         -I"${SCRIPT_DIR}/../../include" \
         -o ./pch.hpp.gch ./pch.hpp
   )
elif [ "$OS_TYPE" = "Linux" ]; then
   # Ubuntu/Linuxの場合
   echo "Linux detected"
   (
      cd "${SCRIPT_DIR}"
      /usr/bin/g++-13 -v -std=gnu++23 -x c++-header -fopenmp -fconcepts -pthread -Ofast -march=native \
         -I"${SCRIPT_DIR}" \
         -I"${SCRIPT_DIR}/../include" \
         -I"${SCRIPT_DIR}/../../include" \
         -o ./pch.hpp.gch ./pch.hpp
   )
else
   echo "Unsupported OS: $OS_TYPE"
fi
