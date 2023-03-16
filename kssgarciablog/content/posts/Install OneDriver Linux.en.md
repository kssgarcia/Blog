---
title: "Install OneDriver Linux"
date: 2023-02-10T11:48:52-05:00
draft: false
---

This post is an installation's guide for OneDriver on Ubuntu or Debian distributions.

Almost my entire life I have used a computer with Windows operating system, some mouth ago I decide to try Linux operating system, just for curiosity, so I made a dual boot in my laptop giving me the possibility to use Linux and Windows, I decided to make this because there are some programs which aren't available on Linux. After some research and post read about Linux and how this operating system is better for programming I exported all my work flow from Windows to Linux, I used Ubuntu distribution for facility. After some weeks I had a problem when I changed to other computers, which in most they has Windows, this was a difficulty for me because I didn't have some tool that allowed me to connect easily the two different system. So I started to search some cloud tool that let me do it, and I found a perfect tool call [OneDriver](https://github.com/jstaf/onedriver#onedriver) easy to use and install, in specific this program is a native Linux file system for Microsoft OneDrive.

Installation:

``` 
$ sudo add-apt-repository ppa:jstaf/onedriver
$ sudo apt update
$ sudo apt install onedriver
```

After make the execution of these commands it will add a `.desktop` file in `/usr/share/applications/` folder. After this, open the [OneDriver](https://github.com/jstaf/onedriver#onedriver) app, setting up the folder where your OneDrive files will remain and make sync with your Microsoft OneDrive account.

Using this program will be easier to keep your work flow when it is necessary to change to other operating systems or computer.
