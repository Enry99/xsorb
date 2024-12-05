xsorb_path=$(pwd)

if !(grep -q ${xsorb_path} ~/.bashrc)
then
echo "export PATH=${xsorb_path}/xsorb/bin:\${PATH}" >> ~/.bashrc
echo "export PYTHONPATH=${xsorb_path}:\${PYTHONPATH}" >> ~/.bashrc
else
echo "xsorb already found on PATH. Please remove it from .bashrc if you want to change the script location."
fi

chmod +x xsorb/bin/xsorb
chmod +x xsorb/bin/xsorb-ml-opt