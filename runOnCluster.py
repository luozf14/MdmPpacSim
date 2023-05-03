from pyspark import StorageLevel, SparkFiles, SparkContext, SparkConf
import subprocess,time,json

def run_sim(index) :
  with open("config.json") as infile :
    config = json.load(infile)
  config["ProcessNumber"] = index
  with open("config_{0}.json".format(index),"w") as outfile :
    json.dump(config,outfile,indent = 4)
  subprocess.call(["./MdmPpacSim","config_{0}.json".format(index)])
  subprocess.call(["mv","Stage1~{0}.root".format(index),"/hdfs/user/luozf/MdmPpacSimResults/"])
#   subprocess.call(["/usr/local/hadoop/bin/hdfs","dfs","-chown","rogachev","/user/rogachev/TexAtSimResults/alpha/fe/MM_{0}.root".format(index)])

if __name__ == "__main__" :
  sconf = SparkConf().setAppName("MdmPpacSim")
  sc = SparkContext(conf=sconf)
    
  sc.addFile("build/MdmPpacSim")
  sc.addFile("run.mac")
  sc.addFile("config/config.json")

  distData = sc.parallelize(range(0,1),1)

  distData.foreach(lambda x: run_sim(x))
