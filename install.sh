xsorb_path=`pwd`

if !(grep -q ${xsorb_path} ~/.bashrc) 
then
echo "export PATH=${xorb_path}:\${PATH}" >> ~/.bashrc
else
echo "xsorb already found on PATH. Please remove it from .bashrc if you want to change the script location."
fi