---
title: "Instalación de OneDriver en Linux"
date: 2023-02-10T11:48:52-05:00
draft: false
---

Esta publicación es una guía de instalación para OneDriver en distribuciones Ubuntu o Debian.

Casi toda mi vida he usado una computadora con sistema operativo Windows, hace algún tiempo decidí probar el sistema operativo Linux, solo por curiosidad, así que hice un arranque dual en mi computadora portátil dándome la posibilidad de usar Linux y Windows, decidí hacer esto porque hay algunos programas que no están disponibles en Linux. Después de investigar un poco y leer sobre Linux y cómo este sistema operativo es mejor para la programación, exporté todo mi flujo de trabajo de Windows a Linux, utilicé la distribución de Ubuntu para la instalación. Después de algunas semanas tuve un problema cuando cambiaba a otras computadoras, que en la mayoría tienen Windows, esto fue una dificultad para mí porque no tenía alguna herramienta que me permitiera conectar fácilmente los dos sistemas diferentes. Así que comencé a buscar alguna herramienta en la nube que me permitiera hacerlo, y encontré una herramienta perfecta llamada OneDriver fácil de usar e instalar, en concreto este programa es un sistema de archivos nativo de Linux para Microsoft OneDrive.

Instalación:

``` 
$ sudo add-apt-repository ppa:jstaf/onedriver
$ sudo apt update
$ sudo apt install onedriver
```

Después de realizar la ejecución de estos comandos se agregará un archivo `.desktop` en la carpeta `/usr/share/applications/`. Posteriormente de esto, abra la aplicación [OneDriver](https://github.com/jstaf/onedriver#onedriver), configure la carpeta donde permanecerán sus archivos de OneDrive y realice la sincronización con su cuenta de Microsoft OneDrive.


El uso de este programa será más fácil para mantener su flujo de trabajo cuando sea necesario cambiar a otros sistemas operativos o computadora.
