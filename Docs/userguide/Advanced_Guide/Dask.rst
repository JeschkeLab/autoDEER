Using with Dask
=====================

Processing data with DeerLab can take a considerable period of time. When this
is happening the Juptyer Notebook is "busy" and unable to perform other task,
so it is advisable to place run computation on a seperate worker thread. 

For large scale problems, this can even be used to run computation on remote 
computers. 

Instalation
----------------------


Basic Usage
----------------------

The first step is to create a client. By default, this is done locally::
    
    from dask.distributed import Client
    client = Client()
    print(client)

It is simplets to create a wrapper function and then to submit said wrapper
to Dask.::

    def test():
        _, _, _, fit = aD.std_deer_analysis(
            BenchB_5pDEER_2hrs.axes,Vre,2.5,2.5,0.2,0.08,100,
            compactness=True,precision="Speed",pathways=[1,2,3,4,5,6,7,8],mask=mask)
        return fit
    
    test_job = client.submit(test)

This then returns a future, that can be accessed.::
    
    test_job.result()