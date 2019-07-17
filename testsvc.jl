using gRPC:run,gRPCServer,gRPCClient,gRPCClient,close,stub
using Rebugger
include("jlout/testsvc_pb.jl")

# implementations of our test services
function Add(req::BinaryOpReq)
    resp = BinaryOpResp()
    resp.result = req.i1 + req.i2
    @info(req.i1 + req.i2)
    resp
end

function Mul(req::BinaryOpReq)
    resp = BinaryOpResp()
    resp.result = req.i1 * req.i
    resp
end
mutable struct TestRpcController <: ProtoRpcController
    debug::Bool
end
# Utility methods for running the server and client
services = (TestMath(Main),)
srvr = gRPCServer(services,19999)
@async run(srvr)
sleep(1)
# myclient = gRPCClient(19999)
# mystub = stub(myclient,TestMathStub)
# controller = TestRpcController(true)
# write_request(myclient.channel,)
# resp = Add(mystub,controller,BinaryOpReq(i1=100,i2=1),(out)->(out))
# show(resp.result)
# close(srvr.sock)
