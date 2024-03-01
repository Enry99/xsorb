xsorb_path=$(pwd)

if !(grep -q ${xsorb_path} ~/.bashrc)
then
echo "export PATH=${xsorb_path}/xsorbed:\${PATH}" >> ~/.bashrc
echo "export PYTHONPATH=${xsorb_path}:\${PYTHONPATH}" >> ~/.bashrc
else
echo "xsorb already found on PATH. Please remove it from .bashrc if you want to change the script location."
fi

chmod +x xsorbed/xsorb
chmod +x xsorbed/xfrag


source ~/.bashrc
