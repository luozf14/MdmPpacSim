from pyspark import StorageLevel, SparkFiles, SparkContext, SparkConf
import subprocess,time,json

def run_sim(index) :
  # Stage1
  with open("config.json") as infile :
    config = json.load(infile)
  config["ProcessNumber"] = index
  with open("config_stage1~{0}.json".format(index),"w") as outfile :
    json.dump(config,outfile,indent = 4)
  subprocess.call(["./MdmPpacSim","config_stage1~{0}.json".format(index)])

  # Stage2
  config["IsTargetChamber"] = "false"
  with open("config_stage2~{0}.json".format(index),"w") as outfile :
    json.dump(config,outfile,indent = 4)
  subprocess.call(["./MdmPpacSim","config_stage2~{0}.json".format(index)])

#   subprocess.call(["hdfs", "dfs", "-moveFromLocal", "-f", "Stage1~{0}.root".format(index), "/user/luozf/MdmPpacSimResults"])

if __name__ == "__main__" :
  sconf = SparkConf().setAppName("MdmPpacSim")
  sc = SparkContext(conf=sconf)
    
  sc.addFile("build/MdmPpacSim")
  sc.addFile("run.mac")
  sc.addFile("config/config.json")

  distData = sc.parallelize(range(0,1),1)

  distData.foreach(lambda x: run_sim(x))
