xorb_path=`pwd`

cp ~/.bashrc ~/.bashrc.bk

if !(grep -q ${xorb_path} ~/.bashrc) 
then
echo "export PATH=${xorb_path}:\${PATH}" >> ~/.bashrc
else
echo "xorb already found on PATH. Please remove it from .bashrc if you want to change the script location."
fi