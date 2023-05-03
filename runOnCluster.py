from pyspark import StorageLevel, SparkFiles, SparkContext, SparkConf
import subprocess,time,json

def run_sim(index) :
  with open("config.json") as infile :
    config = json.load(infile)
  config["ProcessNumber"] = index
  with open("config_{0}.json".format(index),"w") as outfile :
    json.dump(config,outfile,indent = 4)
  subprocess.call(["./MdmPpacSim","config_{0}.json".format(index)])
#   subprocess.call(["mv","MM_{0}.root".format(index),"/hdfs/user/rogachev/TexAtSimResults/alpha/fe"])
#   subprocess.call(["/usr/local/hadoop/bin/hdfs","dfs","-chown","rogachev","/user/rogachev/TexAtSimResults/alpha/fe/MM_{0}.root".format(index)])

if __name__ == "__main__" :
  sconf = SparkConf().setAppName("MdmPpacSim")
  sc = SparkContext(conf=sconf)
    
  sc.addFile("MdmPpacSim")
  sc.addFile("config/config.json")
  sc.addFile("macros/run.mac")

  distData = sc.parallelize(range(0,5),5)

  distData.foreach(lambda x: run_sim(x))