from pyspark import StorageLevel, SparkFiles, SparkContext, SparkConf
import subprocess,time,json,os

def run_sim(index) :
  # Stage1
  with open("config.json") as infile :
    config = json.load(infile)
  config["ProcessNumber"] = index
  config["IsTargetChamber"] = True
  with open("config_stage1~{0}.json".format(index),"w") as outfile :
    json.dump(config,outfile,indent = 4)
  subprocess.call(["source /usr/local/root6/bin/thisroot.sh;", "./MdmPpacSim config_stage1~{0}.json".format(index)],shell=True)

  # Stage2
  config["IsTargetChamber"] = False
  with open("config_stage2~{0}.json".format(index),"w") as outfile :
    json.dump(config,outfile,indent = 4)
  subprocess.call(["source /usr/local/root6/bin/thisroot.sh;", "./MdmPpacSim config_stage2~{0}.json".format(index)],shell=True)

  subprocess.call(["rm", "-rf", "config_stage1~{0}.json".format(index)],shell=True)
  subprocess.call(["rm", "-rf", "config_stage2~{0}.json".format(index)],shell=True)

if __name__ == "__main__" :
  sconf = SparkConf().setAppName("MdmPpacSim")
  sc = SparkContext(conf=sconf)

  # Below are only for run on cluster, should be commented when run on local
  current_path = os.getcwd() 
  subprocess.call('cp build/MdmPpacSim .',shell=True)
  subprocess.call('cp config/config.json .',shell=True)


  distData = sc.parallelize(range(0,2),2)

  distData.foreach(lambda x: run_sim(x))
  subprocess.call('rm MdmPpacSim',shell=True)
  subprocess.call('rm config.json',shell=True)
