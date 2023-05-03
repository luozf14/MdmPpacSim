from pyspark import StorageLevel, SparkFiles, SparkContext, SparkConf
import subprocess,time,json,os

def run_sim(index) :
  subprocess.call("source /usr/local/root6/bin/thisroot.sh",shell=True)
  subprocess.call("source /usr/local/geant4_10.6/bin/geant4.sh",shell=True)
  subprocess.call("source /usr/local/geant4_10.6/share/Geant4-10.6.2/geant4make/geant4make.sh",shell=True)
  
  # Stage1
  with open("config.json") as infile :
    config = json.load(infile)
  config["ProcessNumber"] = index
  with open("config_stage1~{0}.json".format(index),"w") as outfile :
    json.dump(config,outfile,indent = 4)
  subprocess.call(["./MdmPpacSim","config_stage1~{0}.json".format(index)])

  # Stage2
  config["IsTargetChamber"] = False
  with open("config_stage2~{0}.json".format(index),"w") as outfile :
    json.dump(config,outfile,indent = 4)
  subprocess.call(["./MdmPpacSim","config_stage2~{0}.json".format(index)])

  subprocess.call(["hdfs", "dfs", "-moveFromLocal", "-f", "Stage1~{0}.root".format(index), "/user/luozf/MdmPpacSimResults/Stage1~{0}.root".format(index)])
  subprocess.call(["rm", "-rf", "Stage1~{0}.root".format(index)])

if __name__ == "__main__" :
  sconf = SparkConf().setAppName("MdmPpacSim")
  sc = SparkContext(conf=sconf)

  # Below are only for run on cluster, should be commented when run on local
  current_path = os.getcwd() 
  sc.addFile(current_path+'/build/MdmPpacSim')
  sc.addFile(current_path+'/run.mac')
  sc.addFile(current_path+'/config/config.json')

  distData = sc.parallelize(range(0,1),1)

  distData.foreach(lambda x: run_sim(x))
