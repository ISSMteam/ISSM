def QueueRequirements(queudict, queue, np, time):
    #QUEUEREQUIREMENTS - queue requirements in time, number of cpus, by name of queue.
    #
    #   Usage:
    #      QueueRequirements(available_queues, queue_requirements_time, queue_requirements_np, np, time)

    #Ok, go through requirements for current queue:
    try:
        rtime = queudict[queue][0]
    except KeyError:
        raise Exception('QueueRequirements error message: availables queues are ' + queudict.keys)

    if time <= 0:
        raise Exception('QueueRequirements: time should be a positive number')
    if time > rtime:
        raise Exception('QueueRequirements: time should be < ' + str(rtime) + ' for queue: ' + queue)

    #check on np requirements
    if np <= 0:
        raise Exception('QueueRequirements: np should be a positive number')
