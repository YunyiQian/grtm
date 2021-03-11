module comvar
  save
  integer  :: num(10,4), flag(10,4)
  integer  :: flagrp(10), flagip(10)
  real*8   :: rp(10,3), ip(10,3)
  real*8   :: pre(10,5,2), tre(10,5,2), pim(10,5,2),tim(10,5,2)
  real*8   :: sumreal(10), sumimag(10)
end module comvar  