import os
from subprocess import call

if __name__ == '__main__':

    tmp =[]
    for i in os.listdir("/home/ubuntu"):
        if i.endswith("_L1") or i.endswith("_L2"):
            for j in os.listdir("/home/ubuntu/{}".format(i)):
                if j.endswith(".genes.results"):
                    tmp.append("/home/ubuntu/{}/{}".format(i,j))
    print(tmp)
            # if os.path.exists("/home/ubuntu/{}/RSEM.genes.results".format(i)):
            #     marker= "{}-{}-{}".format(i.strip().split("_")[0],i.strip().split("_")[1],i.strip().split("_")[2])
            #     call("mv /home/ubuntu/{}/RSEM.genes.results /home/ubuntu/{}/{}.genes.results".format(i,i,marker),shell=True)

